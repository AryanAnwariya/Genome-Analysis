# GM12878 Compartment Knowledge Graph

## Overview

This component organises genomic bins and their
compartment classifications in a SQLite database.

It provides a way to look up genomic regions and
compare their classifications from different methods.

## Data

- Cell line: GM12878
- Resolution: 100 kb
- Initial database: 28,823 genomic bins
- Initial build reported 373,478 classification records
- Sources include published compartment annotations
  and graph embedding clustering results.

The actual database tables and record counts are
documented in validation.txt.

## Implementation

compartment_kg.py contains the database implementation.

## Results

- validation.txt: database schema, record counts,
  and integrity verification
- sample_query.txt: genomic region lookup example

## Relationship to the Hi-C Query Processor

The compartment knowledge graph organises annotations
for genomic regions.

The hictile query processor retrieves contact values
directly from Hi-C matrices.

These components are separate in the current prototype.
Future work may integrate them.

## Limitations

The initial knowledge graph is a structured compartment
annotation database, not a completed biological knowledge
graph covering all genomic relationships.

The SQLite database and large input datasets remain on HPC.
