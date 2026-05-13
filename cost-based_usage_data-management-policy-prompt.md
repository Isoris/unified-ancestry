# Cost-Based Data Management Policy (paste-into-another-chat prompt)

**What this file is.** A portable instruction set for an AI assistant
(Claude / similar) to audit a research codebase and classify every
computed result by **cost × reuse**, then route each to the right
storage tier. Paste this whole file into a fresh chat with the
instruction *"audit this repo per the policy below"*.

---

## The one-line rule

> Cheap results are streamed live. Expensive reusable results are
> cached by hash. Very expensive shared layers are precomputed.

## Latency thresholds (default classifier)

| compute time         | tier                                    |
|----------------------|-----------------------------------------|
| < 50 ms              | live compute                            |
| 50 – 500 ms          | live + memory cache                     |
| 0.5 – 5 sec          | tmp / project cache if repeated         |
| > 5 sec              | project cache or precompute             |
| minutes – hours      | always precompute / shared layer        |

## Decision matrix (data class → tier)

| data class                                       | tier                |
|--------------------------------------------------|---------------------|
| UI formatting / colors / labels                  | live                |
| simple filter / sort / rank                      | live                |
| tooltip from already-loaded JSON                 | browser memory only |
| small summary from one result.json               | live                |
| FST for one small interval                       | live (maybe)        |
| FST for 50 candidates                            | disk/project cache  |
| local PCA coordinates                            | precompute          |
| ancestry Q windows                               | precompute          |
| ROH genome-wide summaries                        | precompute          |
| π / dXY / Tajima's D window tracks               | project cache       |
| random parameter exploration                     | tmp cache or none   |
| final reviewed / manuscript result               | committed cache     |

## The five storage tiers

1. **Live compute** — no storage. Recompute every request.
2. **Memory / session cache** — browser or server RAM, ephemeral.
3. **Disk / project cache** — hashed, e.g.
   `popstats_cache/computed/<analysis_type>/<sha256_or_fnv64>/result.json`.
   Cache key = canonicalized request + input-file mtimes.
4. **Precompute** — shared layer, written once into `popstats_cache/`.
   Atlas sessions read; never re-compute interactively.
5. **Committed / saved** — git or release package. Reviewed, signed off.

## Mapping by usage pattern (cross-axis)

| reuse pattern                              | typical tier           |
|--------------------------------------------|------------------------|
| computed once, used everywhere             | precompute             |
| computed often, deterministic              | disk/project cache     |
| computed often, depends on session params  | memory/session cache   |
| computed rarely, exploratory               | tmp cache or no cache  |
| computed in the moment, displayed only     | live                   |

## Audit procedure (what the assistant should do)

For each distinct computed-result type in the target repo:

1. **Find** — grep for compute-emitting code paths (functions that
   write JSON / TSV / RDS / pickle; HTTP endpoints; scheduled jobs).
2. **Estimate cost** — read the code, identify input size + algorithm
   complexity. If unsure, run on representative inputs and time it.
3. **Estimate reuse** — does the same result get requested again by
   the same user in the same session? Across sessions? Across users?
4. **Classify** — pick a tier from the latency table OR the decision
   matrix; cross-check with the reuse pattern.
5. **Recommend a storage path** — relative path under
   `popstats_cache/`, `<repo>/cache/`, etc.
6. **Specify an invalidation rule** — input-file mtimes, content hash,
   or "manual / never".
7. **Note current behavior** — is the code already doing this? If
   not, what's the mismatch?

## Output the audit as

`cost_audit.tsv` — one row per result type, columns:

| col | meaning |
|---|---|
| `result_type`    | short id (e.g. `fst_per_window`, `local_Q`) |
| `description`    | one-sentence what-it-is |
| `where_emitted`  | file:line of the compute / writer |
| `cost_estimate`  | rough wall-time on representative input |
| `reuse_pattern`  | once / session / cross-session / shared |
| `recommended_tier` | live / memory / disk / precompute / committed |
| `recommended_path` | concrete path |
| `invalidation`   | mtime / content-hash / manual / never |
| `current_behavior` | what the code does today |
| `mismatch`       | "ok" if current_behavior matches recommended; else describe |
| `priority`       | low / medium / high to fix |

Plus `cost_audit.md` — short prose summary:
- which tiers are over-/under-used,
- top-3 highest-priority mismatches,
- proposed migration order.

## Apply-it cheatsheet (paste into Claude as the task)

> Audit this repo against `cost-based_usage_data-management-policy-prompt.md`.
> Produce `cost_audit.tsv` and `cost_audit.md` at repo root. Cite
> file:line for every entry. Do not modify any compute code; the audit
> is doc-only on this pass.

## Why this matters (one paragraph for new readers)

A research codebase tends to drift toward "cache everything, precompute
everything" because each new result feels expensive in isolation. That
inflates storage, slows builds, and ages out fast when inputs change.
The opposite extreme — "compute everything live" — kills interactive
latency for any non-trivial scan. The policy above splits the
difference by **cost × reuse**: cheap results stay live (no storage
debt), expensive deterministic results get hashed disk cache
(invalidate-on-input-change for free), and the genuinely heavy shared
layers (ancestry Q, PCA, ROH, callable masks) get precomputed once and
read everywhere.

The atlas streams live where possible, hits the cache where it must,
and treats precompute as the load-bearing layer behind it.
