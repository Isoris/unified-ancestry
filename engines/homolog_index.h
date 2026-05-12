// =============================================================================
// homolog_index.h — shared binary record definitions for the homolog index.
//
// Layout (little-endian, no padding inside records — packed):
//   [0]                          IndexHeader (fixed-size)
//   [off_gene_index]             GeneRecord[n_genes]   — sorted alphabetically by gene name
//   [off_hits]                   HitRecord[n_hits]     — clustered by gene, then bitscore desc
//   [off_string_pool]            char[string_pool_size] — NUL-terminated strings,
//                                                          referenced by *_off fields
//
// Schema documented in engines/schemas/homolog_index.binary.schema.json.
// =============================================================================

#ifndef HOMOLOG_INDEX_H
#define HOMOLOG_INDEX_H

#include <stdint.h>

// "HOLINDX\x01" little-endian = 0x0158444E494C4F48
#define HOLINDX_MAGIC   0x0158444E494C4F48ULL
#define HOLINDX_VERSION 1

#define HOLINDX_SRC_DIAMOND  0
#define HOLINDX_SRC_MINIPROT 1

typedef struct __attribute__((packed)) {
    uint64_t magic;
    uint32_t version;
    uint32_t n_genes;
    uint64_t n_hits;
    uint64_t string_pool_size;
    uint64_t header_size;
    uint32_t gene_record_size;
    uint32_t hit_record_size;
    uint64_t off_gene_index;
    uint64_t off_hits;
    uint64_t off_string_pool;
    uint8_t  reserved[16];
} IndexHeader;

typedef struct __attribute__((packed)) {
    uint32_t name_off;     // offset into string pool
    uint32_t name_len;     // length excluding NUL (informational)
    uint64_t hits_off;     // index into hit array (NOT a byte offset)
    uint32_t n_hits;
    uint32_t reserved;
} GeneRecord;

typedef struct __attribute__((packed)) {
    uint32_t target_off;   // string pool offset of target id (subject for DIAMOND, contig for miniprot)
    uint32_t source_off;   // string pool offset of source label (e.g. species name)
    int32_t  qstart, qend; // 1-based, inclusive (as DIAMOND emits) or 0-based half-open for PAF — preserved as-is
    int32_t  tstart, tend;
    float    identity;     // fractional 0..1 (DIAMOND pident/100; miniprot matches/alnlen)
    float    bitscore;     // -1 if not available
    double   evalue;       // -1 if not available
    uint32_t alnlen;
    uint16_t n_segments;   // 1 for DIAMOND; could be parsed from miniprot cs:Z: tag later
    uint8_t  source_type;  // HOLINDX_SRC_*
    uint8_t  strand;       // '+' / '-' / 0
} HitRecord;

#endif
