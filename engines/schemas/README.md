# engines/schemas — JSON contracts for the popstats local server

Two-half contract per statistic:

- **Activator schema** = `<stat>.activator.schema.json`. The JSON request the
  client posts to the popstats server to invoke a compute. Names the statistic,
  inputs (paths/refs), scope (chrom / windows), and parameters.
- **Extractor schema** = `<stat>.extractor.schema.json`. The JSON response the
  server returns. Carries `schema_version`, `metadata`, and the per-window /
  per-region records.

The TSVs the C binaries currently emit are an intermediate that the dispatcher
converts to the extractor JSON shape. The C binaries are statistic-agnostic
plumbing; the schemas pin down the wire format.

## Per-statistic schemas in this dir

| Statistic | Engine | Activator | Extractor |
|---|---|---|---|
| codon_stats | engines/codon_stats | (TSV-only legacy) | codon_stats.output.schema.json |
| homolog_index lookup | engines/homolog_index_query | — (CLI args) | homolog_index.query.output.schema.json |
| homolog atlas HTTP | engines/homolog_atlas_server | (URL query params) | homolog_atlas_server.api.schema.json |
| xpehh | engines/xpehh | xpehh.activator.schema.json | xpehh.extractor.schema.json |
| iHS | engines/iHS | iHS.activator.schema.json | iHS.extractor.schema.json |
| outlier_scan | engines/outlier_scan | outlier_scan.activator.schema.json | outlier_scan.extractor.schema.json |
| candidate_vs_flanks | engines/candidate_vs_flanks | candidate_vs_flanks.activator.schema.json | candidate_vs_flanks.extractor.schema.json |
| (umbrella) | popstats server | popstats_server.activator.schema.json | popstats_server.extractor.schema.json |

The umbrella schema is a `oneOf` over each per-statistic activator/extractor;
clients can target either the umbrella (with `statistic` field as the
discriminator) or the per-statistic schemas directly.

## Status notes

This dir captures the activator/extractor contracts for the NEW engines I'm
adding (xpehh, iHS, +). Older engines (codon_stats, region_popstats, etc.)
still rely on TSV outputs only; their activator schemas are planned but not
yet written — the dispatcher can validate activator inputs by hand for those
until they catch up. Filling in those gaps is a separate ticket.

## Caching is the server's job

`engines/popstats_dispatch` is stateless: activator JSON in → engine run →
extractor JSON streamed out. There is no built-in cache.

Whatever HTTP server fronts the dispatcher is the right layer to memoize
responses, because the server already:

- holds the raw activator JSON before it's parsed (cheap content hash);
- knows the deployment-specific storage backend (Redis, filesystem, in-memory
  LRU, …);
- can enforce per-user or per-route TTL / eviction policy;
- can decide cache-key canonicalization (drop `request_id`, sort keys, etc.).

Caching strategy recommended at the server layer:

1. Compute `key = hash(canonicalize(activator - {request_id}))`.
2. Resolve referenced input file mtimes (`inputs.beagle`, `inputs.windows_tsv`,
   `inputs.candidates_bed`, …); include those mtimes in the cache key OR
   invalidate the entry when any input is newer than the stored response.
3. Cache the extractor JSON verbatim. It's already the wire-format payload.
4. On hit: stream the cached bytes to the client; on miss: pipe the activator
   into `popstats_dispatch`, capture stdout, store, stream.

If at some future point this turns out to be the wrong split (e.g. CLI users
also want caching), `--cache_dir` can be added to the dispatcher later as
opt-in without breaking the stateless default.
