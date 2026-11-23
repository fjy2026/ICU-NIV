"""MedCortex 运行追踪的检查点（由 AnalysisTrace 生成，请勿修改）。

用法：from mc_trace import mc_checkpoint
      df = mc_checkpoint(df, "排除年龄 < 18 岁", id="subject_id")
      mc_checkpoint(df, "排除妊娠", id="SEQN", keep=df["in_cohort"])   # 只打标记不删行时：只数标记为真的
      rx = mc_checkpoint(rx, "处方药表", id="SEQN", side=True)          # 单独读入、要并回主表的另一张表
读入数据后打一个，每次筛选 / 合并后再打一个。人数由程序在真实数据上数，写进运行记录；
不在 AnalysisTrace 里运行时什么都不做。支持 pandas 与 polars 数据框。
"""
import hashlib
import json
import os
import sys


def _key(value):
    # 与 R 版 .mc_id_hash 同一套写法：整数值写成整数（-0 也是 0），其余小数按 17 位有效数字，逻辑值写成 TRUE/FALSE；NaN、空白字符串当缺失
    if value is None:
        return None
    if isinstance(value, bool):
        return 'TRUE' if value else 'FALSE'
    if isinstance(value, float):
        if value != value:
            return None
        return str(int(value)) if value.is_integer() else format(value, '.17g')
    text = str(value)
    return text if text.strip() else None


def _true(value):
    try:
        return bool(value) and value == value
    except TypeError:  # pd.NA 之类的缺失
        return False


def _record_env():
    # 脚本进程里实际加载的包与版本：复跑结果不同时据此判断是不是环境变了（每次覆盖，记到最后一个检查点为止）
    path = os.environ.get('MC_TRACE_ENV')
    if not path:
        return
    packages = {}
    for name, module in list(sys.modules.items()):
        version = getattr(module, '__version__', None)
        if '.' not in name and not name.startswith('_') and isinstance(version, str):
            packages[name] = version
    with open(path, 'w', encoding='utf-8') as stream:
        json.dump(dict(sorted(packages.items())[:80]), stream, ensure_ascii=False)


def mc_checkpoint(df, label, id=None, keep=None, side=False):
    path = os.environ.get('MC_TRACE_CHECKPOINTS')
    if not path:
        return df
    _record_env()
    view = df
    if keep is not None:
        mask = [_true(k) for k in (keep.to_list() if hasattr(keep, 'to_list') else list(keep))]
        if hasattr(df, 'isna'):
            view = df[mask]
        else:
            import polars
            view = df.filter(polars.Series(mask))
    record = dict(label=str(label), rows=int(len(view)), columns=len(list(view.columns)))
    if side:
        record['side'] = True
    if id is not None:
        values = view[id]
        values = values.dropna() if hasattr(values, 'dropna') else values.drop_nulls()
        keys = sorted({k for k in (_key(v) for v in values.to_list()) if k is not None})
        record.update(id=str(id), n_id=len(keys),
                      id_hash=hashlib.md5(''.join(k + '\n' for k in keys).encode('utf-8')).hexdigest() if keys else None)
        if keys:
            # 保留的编号（只存编号）：核验时可从被排除的人里抽查，排除的人 = 上一检查点编号 - 本检查点编号
            folder = os.path.join(os.path.dirname(path), 'ids')
            os.makedirs(folder, exist_ok=True)
            target = os.path.join(folder, record['id_hash'] + '.txt.gz')
            if not os.path.exists(target):
                import gzip
                with gzip.open(target, 'wt', encoding='utf-8', newline='') as stream:  # 不让 Windows 把 \n 写成 \r\n：与 R 侧及文件名里的哈希一致
                    stream.write(''.join(k + '\n' for k in keys))
    nulls = view.isna().sum().to_dict() if hasattr(view, 'isna') else view.null_count().row(0, named=True)
    record['missing'] = {str(k): int(v) for k, v in list(nulls.items())[:500]}
    # 分类列（不超过 10 种取值）各取值的个数：交付说明里写的「男性 N 人」「死亡 N 例」由它核对；表太大时不数
    if len(view) <= 2000000:
        levels = {}
        for column in list(view.columns)[:500]:
            if column == id:
                continue
            series = view[column]
            series = series.dropna() if hasattr(series, 'dropna') else series.drop_nulls()
            pandas = hasattr(series, 'nunique')
            if (series.nunique() if pandas else series.n_unique()) > 10:
                continue
            counts = {}
            for value, count in (series.value_counts().items() if pandas else series.value_counts().iter_rows()):
                key = _key(value)
                if key is not None and count:
                    counts[key] = counts.get(key, 0) + int(count)
            if counts:
                levels[str(column)] = counts
        record['levels'] = levels
        # 数值列（超过 10 种取值）的最小、5% 分位、95% 分位、最大：收尾「取值体检」据此点名 9999 这类编码值、单位错
        ranges = {}
        for column in list(view.columns)[:500]:
            if column == id or str(column) in levels:
                continue
            series = view[column]
            pandas = hasattr(series, 'nunique')
            if not (series.dtype.kind in 'iuf' if pandas else series.dtype.is_numeric()):
                continue
            if pandas:
                series = series.dropna().astype(float)
                series = series[(series != float('inf')) & (series != float('-inf'))]
                if len(series):
                    p5, p95 = series.quantile([0.05, 0.95]).to_list()
            else:
                import polars
                series = series.drop_nulls().cast(polars.Float64)
                series = series.filter(series.is_finite())
                if len(series):
                    p5, p95 = series.quantile(0.05, 'linear'), series.quantile(0.95, 'linear')
            if len(series):
                ranges[str(column)] = [float(series.min()), float(p5), float(p95), float(series.max())]
        record['ranges'] = ranges
    frame = sys._getframe(1)
    record.update(file=os.path.abspath(frame.f_code.co_filename), line=frame.f_lineno)
    # 在函数库里调用的：顺调用栈找到主脚本里那一行，追溯面板据此高亮主脚本
    main, outer = os.path.abspath(sys.argv[0]) if sys.argv and sys.argv[0] else '', frame
    while outer is not None and os.path.abspath(outer.f_code.co_filename) != main:
        outer = outer.f_back
    if outer is not None and outer is not frame:
        record['entry_line'] = outer.f_lineno
    with open(path, 'a', encoding='utf-8') as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + '\n')
    return df
