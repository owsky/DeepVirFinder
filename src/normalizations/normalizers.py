from normalizations import count, divmax, log, min_max, z_score, raw_count, z_score_count

normalizers = {
    "count": count.normalize_sequence,
    "divmax": divmax.normalize_sequence,
    "log": log.normalize_sequence,
    "min_max": min_max.normalize_sequence,
    "z_score": z_score.normalize_sequence,
    "raw_count": raw_count.normalize_sequence,
    "z_score_count": z_score_count.normalize_sequence
}