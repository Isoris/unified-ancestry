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
