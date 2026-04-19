"""Shared condition metadata for analyze scripts (aligned with generate_rank_lists.py)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator


@dataclass(frozen=True)
class ConditionSpec:
    """One experimental arm; `condition_id` matches the .rnk filename stem."""

    condition_id: str
    log2fc_col: str
    neglogp_col: str


# Order matches output/*.rnk naming from generate_rank_lists.py
CONDITIONS: tuple[ConditionSpec, ...] = (
    ConditionSpec("αCD3_5min", "Log2FC_αCD3_5min", "-Log_p-value_αCD3_5min"),
    ConditionSpec("αCD3_αCD226_5min", "Log2FC_αCD3_αCD226_5min", "-Log_p-value_αCD3_αCD226_5min"),
    ConditionSpec("αCD3_αICOS_5min", "Log2FC_αCD3_αICOS_5min", "-Log_p-value_αCD3_αICOS_5min"),
    ConditionSpec(
        "αCD3_αCD2_5min",
        "Log2FC_αCD3_αCD2_5min)",  # typo in data.tsv header
        "-Log_p-value_αCD3_αCD2_5min",
    ),
    ConditionSpec("αCD3_10min", "Log2FC_αCD3_10min", "-Log_p-value_αCD3_10min"),
    ConditionSpec(
        "αCD3_αCD226_10min",
        "Log2FC_αCD3_αCD226_10min",
        "-Log_p-value_αCD3+αCD226_10min",
    ),
    ConditionSpec("αCD3_αICOS_10min", "Log2FC_αCD3_αICOS_10min", "-Log_p-value_αCD3_αICOS_10min"),
    ConditionSpec(
        "αCD3_αCD2_10min",
        "Log2FC_αCD3_αCD2_10min",
        "-Log_p-value_αCD3+αCD2_10min",
    ),
)


@dataclass(frozen=True)
class TemporalArm:
    """5 min vs 10 min within the same co-stimulation background."""

    arm_id: str
    t5_condition_id: str
    t10_condition_id: str


TEMPORAL_ARMS: tuple[TemporalArm, ...] = (
    TemporalArm("aCD3_only", "αCD3_5min", "αCD3_10min"),
    TemporalArm("aCD3_aCD226", "αCD3_αCD226_5min", "αCD3_αCD226_10min"),
    TemporalArm("aCD3_aICOS", "αCD3_αICOS_5min", "αCD3_αICOS_10min"),
    TemporalArm("aCD3_aCD2", "αCD3_αCD2_5min", "αCD3_αCD2_10min"),
)


def condition_by_id(cid: str) -> ConditionSpec:
    for c in CONDITIONS:
        if c.condition_id == cid:
            return c
    raise KeyError(f"Unknown condition_id: {cid!r}")


def iter_condition_ids() -> Iterator[str]:
    yield from (c.condition_id for c in CONDITIONS)
