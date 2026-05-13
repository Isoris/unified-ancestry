# registries/

Single source of truth for what this repo exposes, indexed in
`manifest.json` (the registry-of-registries). Read that first if
you're a downstream consumer; the per-registry files below are
just the data.

## Inventory

| file | what | schema |
|---|---|---|
| `manifest.json`             | umbrella — lists every registry below, with paths + schemas | (self-describing, `registry_manifest_version: 1`) |
| `analysis_capability.json`  | what analyses this repo ships (engine-side fields only) | `engines/schemas/analysis_capability.entry.schema.json` |
| `interval_registry.tsv`     | intervals at all scales (genome / chrom / SNP windows / candidates) | `engines/schemas/interval_registry.schema.json` |
| `sample_subsets.tsv`        | named sample groupings → paths to ID-per-line files | `engines/schemas/sample_subsets.schema.json` |
| `cov_registry.tsv`          | PCAngsd .cov files indexed by chrom/scope | `engines/schemas/cov_registry.schema.json` |
| `build_registries.py`       | builder for the 3 TSV registries | — |

## Two kinds of registry

**Engine-side manifest (JSON, hand-maintained).**
`analysis_capability.json` is committed to the repo and edited by hand
when an engine is added or its capabilities change. It declares
what each analysis does (analysis_type, provides, requires,
supported_execution_modes, result_schema, entrypoint) but NOT how a
specific deployment wires it (default_policy, cost_class, reuse_scope
live in the consumer server's deployment-policy file, not here).

**Discovery TSVs (auto-built, machine-emitted).**
`interval_registry.tsv`, `sample_subsets.tsv`, `cov_registry.tsv` are
rebuilt by `build_registries.py` from the cohort's config + on-disk
artifacts. Re-run the builder when the cohort or assembly changes.

## Building the TSV registries

```bash
python3 registries/build_registries.py \
  --config 00_ancestry_config.sh \
  --outdir registries/ \
  [--candidates path/to/candidates.tsv]
```

Emits `interval_registry.tsv`, `sample_subsets.tsv`, `cov_registry.tsv`.

## Cross-repo: how the atlas / popstats server consumes this

```
unified-ancestry/registries/manifest.json
       │
       ├──► analysis_capability.json
       │       └─ server adds deployment policy on top
       │           (default_policy, cost_class, reuse_scope)
       │           per cost-based_usage_data-management-policy-prompt.md
       │
       ├──► interval_registry.tsv  ──── filter rows ──► picks the interval(s) for an analysis
       ├──► sample_subsets.tsv     ──── lookup id   ──► passes sample_list_path to the engine
       └──► cov_registry.tsv       ──── lookup chrom─► passes cov_path to PCA-using analyses
```

The schemas live in `engines/schemas/` (not duplicated here) so they
can be `$ref`-imported from sister repos without copying.

## Field ownership reminder

For `analysis_capability.json`, only the **engine-side** fields are
filled in this repo:

| field | who fills |
|---|---|
| `analysis_type`, `module`, `provides`, `requires`, `supported_execution_modes`, `result_schema`, `entrypoint`, `engine_binding_ref` | **engine (this repo)** |
| `default_policy`, `allowed_policies`, `cost_class`, `reuse_scope` | **deployment / consumer server (sister repo)** |

The schema (`analysis_capability.entry.schema.json`) marks deployment
fields as optional. A consumer's deployment registry extends each
engine entry with its policy choice — no edits to this repo needed
when those choices change.

## Adding a new analysis

1. Implement the engine binary (or extend an existing one).
2. Add a row to `analysis_capability.json` with engine-side fields.
3. If the analysis is dispatcher-routed (via `engines/popstats_dispatch`),
   also write `<stat>.activator.schema.json` + `<stat>.extractor.schema.json`
   in `engines/schemas/` and append a `oneOf` `$ref` to the umbrella
   `popstats_server.activator.schema.json` / `.extractor.schema.json`.
4. Update `ENGINES.md` so the human-readable inventory stays in sync.
