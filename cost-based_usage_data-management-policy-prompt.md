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

## How to classify a request — DO NOT hardcode by analysis name

> "FST is cache" / "PCA is precompute" is wrong on principle. The same
> analysis can land in any tier depending on the *instance* — input size,
> candidate count, parameter combination, who's asking. The policy is a
> **classifier**, not a list.

For each request that hits the system, route it by **two registered
properties of the request, not by the analysis name**:

1. **`cost_class`** — measured or declared compute cost of *this
   request shape*. Bucketed by the latency table above.
2. **`reuse_pattern`** — how often *this same request* is expected to
   be issued again. One of:
   - `one_shot`               — never re-requested
   - `session_repeat`         — repeated in the same UI session
   - `cross_session_determ`   — same hash, deterministic, reused across sessions
   - `shared_layer`           — every user / page wants it
   - `exploratory`            — random parameter sweep, rarely hits the same key twice
   - `reviewed_artifact`      — saved, signed off, manuscript-bound

Apply the classifier:

```
cost_class       × reuse_pattern               →  tier
─────────────────────────────────────────────────────────────────────────
cheap            × any                         →  live
medium           × one_shot                    →  live
medium           × session_repeat              →  memory/session cache
medium           × cross_session_determ        →  disk/project cache
expensive        × one_shot / exploratory      →  live or tmp cache
expensive        × cross_session_determ        →  disk/project cache (hashed)
expensive        × shared_layer                →  precompute
very_expensive   × any                         →  precompute (always)
any              × reviewed_artifact           →  committed / saved
```

The **same analysis** (e.g. "FST scan") routes differently in different
contexts:

- FST inside a 5 kb candidate, one user, one click → `medium` × `one_shot` → live
- FST across 50 candidates, scripted batch → `expensive` × `cross_session_determ` → disk cache
- FST for the genome-wide background panel every atlas user reads → `very_expensive` × `shared_layer` → precompute

## The five storage tiers

1. **Live compute** — no storage. Recompute every request.
2. **Memory / session cache** — browser or server RAM, ephemeral.
3. **Disk / project cache** — hashed, e.g.
   `popstats_cache/computed/<analysis_type>/<sha256_or_fnv64>/result.json`.
   Cache key = canonicalized request + input-file mtimes.
4. **Precompute** — shared layer, written once into `popstats_cache/`.
   Atlas sessions read; never re-compute interactively.
5. **Committed / saved** — git or release package. Reviewed, signed off.

## Step 0 — the analysis capability registry

Before the audit can classify anything, the system needs each analysis
to be **self-describing**. The existing "loaded modules" registry only
says `loaded: yes` — that's not enough. The server can't decide whether
to run a request live, cache it, or precompute it if the analysis doesn't
tell it:

- what does this analysis *provide*?
- what inputs does it *require*?
- how expensive is it?
- which execution modes does it *support* (live / memory / disk / precompute / committed)?
- what's its *default* mode?
- what result schema does it *emit*?
- what's the engine *entrypoint*?

The missing piece is an **analysis capability registry** — a small JSON
or TSV loaded at startup. One entry per analysis_type. Minimal shape:

```json
{
  "analysis_type":     "fst_profile",
  "module":            "popstats",
  "provides":          ["fst_profile_v1"],
  "requires":          ["interval", "groups", "cache"],
  "execution_modes": {
    "live":         true,
    "memory_cache": true,
    "disk_cache":   true,
    "precompute":   false,
    "committed":    true
  },
  "default_policy":    "disk_cache",
  "cost_class":        "medium",
  "reuse_scope":       "candidate",
  "result_schema":     "fst_result_v1",
  "entrypoint":        "popstats_fst"
}
```

With this in place the server can ask:
- *"Who provides `fst_profile_v1`?"* — registry answers.
- *"Can it run live?"* — registry says yes/no.
- *"What's the default policy for this analysis?"* — registry says.
- *"Where should the cached result.json go?"* — derived from
  `default_policy` × `analysis_type` × request hash.

Then the cost-based classifier (above) overrides the default per
request shape: e.g. an `fst_profile` request whose request-shape is
genome-wide (shared layer) gets re-routed from `disk_cache` →
`precompute` even though the analysis's default is `disk_cache`.

JSON Schema for one entry:
`engines/schemas/analysis_capability.entry.schema.json` (in this
repo). Other repos should keep the same schema and a sibling
`<repo>/registry/analysis_capability.json` (or .tsv) listing their
analyses.

---

## What the audit produces — per registered module/analysis

The job isn't to slap every "FST" call into the same tier. The job is
to find every **registered module / analysis / endpoint** in the repo,
read off its declared inputs and typical request shape, and *classify
that module's instances*. A module can have multiple rows in the audit
if it's called in materially different shapes.

For each module, the audit walks its known callers and produces one
row per (module × request-shape) pair.

If the analysis capability registry above exists, the audit reads it
first as the source of truth. If it doesn't, **building it is the
first audit deliverable** — the audit produces the registry as a
side-effect.

## Audit procedure (what the assistant should do)

For each **registered module / analysis / endpoint** in the target
repo, AND for each materially different **request shape** the module
is called with:

1. **Find the module** — registry file (e.g. `registries/*.tsv`,
   `engines/Makefile`, dispatcher routes, HTTP routes). Cite file:line.
2. **Find the request shapes** — what inputs and parameter ranges does
   the module get called with in practice? Read the dispatchers,
   atlas pages, scheduled jobs. Group similar shapes together.
3. **Estimate `cost_class`** — per request shape. Read the code,
   identify input size × algorithm complexity. If unsure, run on
   representative inputs and time it. Bucket per the latency table.
4. **Estimate `reuse_pattern`** — per request shape. One of:
   `one_shot` / `session_repeat` / `cross_session_determ` /
   `shared_layer` / `exploratory` / `reviewed_artifact`.
5. **Classify** — apply the cost × reuse table above. The same module
   can produce several rows in the audit if it's called in different
   shapes that route to different tiers.
6. **Recommend a storage path** — concrete relative path. Examples:
   `popstats_cache/precomputed/<module>/<chrom>/result.json` for
   precompute, `popstats_cache/computed/<module>/<request_hash>/result.json`
   for hashed cache, `null` for live.
7. **Specify an invalidation rule** — input-file mtimes, content
   hash, registry version, or "manual / never".
8. **Note current behavior** — does the code already route this
   request shape to that tier? If not, describe the mismatch.

## Output the audit as

`cost_audit.tsv` — one row per (module × request-shape) pair, columns:

| col | meaning |
|---|---|
| `module`           | registry id of the analysis / endpoint (e.g. `fst_window`, `local_Q`, `pi_NS`) |
| `request_shape`    | what's different about this call vs other rows of the same module (e.g. "single 5 kb interval", "50-candidate batch", "genome-wide background") |
| `description`      | one-sentence what-it-is for this shape |
| `where_registered` | file:line of the registry / route / dispatcher binding |
| `cost_class`       | cheap / medium / expensive / very_expensive |
| `cost_estimate`    | rough wall-time on representative input |
| `reuse_pattern`    | one_shot / session_repeat / cross_session_determ / shared_layer / exploratory / reviewed_artifact |
| `recommended_tier` | live / memory / disk / precompute / committed (from the classifier) |
| `recommended_path` | concrete path or `null` for live |
| `invalidation`     | mtime / content-hash / registry-version / manual / never |
| `current_behavior` | what the code does today for this shape |
| `mismatch`         | "ok" if current_behavior matches recommended; else describe |
| `priority`         | low / medium / high to fix |

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
