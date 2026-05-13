# Cost-Based Cache Policy (paste-into-another-chat prompt)

**What this file is.** A short, portable instruction set for an AI
assistant. Paste this whole file into a fresh chat with the
instruction *"audit this repo per the policy below"*. The aim is
**make the atlas go fast** by caching the right things, not
overengineering.

## The rule

> **Cheap results stream live. Expensive deterministic results get
> cached on disk. Very heavy shared layers get precomputed.**

That's it. Three tiers.

## Latency → tier

| compute time     | tier                  | where it lives                                              |
|------------------|-----------------------|--------------------------------------------------------------|
| < 100 ms         | **live**              | nowhere; recompute every time                                |
| 100 ms – 5 sec   | **disk cache**        | `cache/<analysis>/<request_hash>/result.json`                |
| > 5 sec, shared  | **precompute**        | `cache/precomputed/<analysis>/...` written by offline jobs   |

Default: **disk cache**. Promote to precompute when many users hit
the same result. Demote to live when caching costs more than recompute.

## Cache key

`request_hash = fnv64(canonicalize(request_json - {request_id}))`

Invalidate when any input file's mtime is newer than the cache file's
mtime. That's the whole invalidation rule.

## Self-describing analyses

Each toolkit repo ships a manifest at `registries/analysis_capability.json`
listing every analysis it exposes. Each entry says: `analysis_type`,
`module`, `provides`, `requires`, `supported_execution_modes`,
`result_schema`, `entrypoint`. Schema:
`engines/schemas/analysis_capability.entry.schema.json`.

The downstream server reads the manifest, picks the tier per the
latency table above, writes results to `cache/<analysis>/<hash>/`.
No further metadata needed.

## What the audit produces

`cost_audit.tsv`, one row per analysis × request-shape:

| col | meaning |
|---|---|
| `analysis_type`    | from the manifest |
| `request_shape`    | what makes this row distinct (small interval / full chrom / genome-wide) |
| `wall_time_estimate` | rough seconds on representative input |
| `recommended_tier` | live / disk_cache / precompute |
| `cache_path`       | concrete path, or `null` for live |
| `current_behavior` | what the code does today |
| `mismatch`         | "ok" if matches; else describe |

Plus `cost_audit.md` — 3 bullets: tier distribution, top mismatch,
one-line migration recommendation.

## Paste this into Claude per repo

> Audit this repo against `cost-based_usage_data-management-policy-prompt.md`.
> Produce `cost_audit.tsv` and `cost_audit.md` at repo root. Doc only;
> don't change compute code.

## One-sentence summary

Use disk cache by default, time things to decide, and don't cache what's
already fast.
