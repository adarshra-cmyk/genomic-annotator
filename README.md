# Genomic Annotator 


## Quick Start

### Installation

```bash
git clone https://github.com/adarshra-cmyk/genomic-annotator.git
cd genomic-annotator
pip install -r requirements.txt
pip install -e .
```

### Basic Usage

```python
from genomic_annotator import AnnotationPipeline

# Initialize pipeline
pipeline = AnnotationPipeline()

# Annotate a single variant
pipeline.run_single_variant("rs238242")

# Process multiple variants
variants = ["rs238242", "rs7527068", "rs142513484"]
results_df = pipeline.run_batch_variants(variants)
print(results_df)

# Process from CSV file
results = pipeline.run_from_csv("my_variants.csv", output_file="results.csv")
```

## Examples

Run the examples to see the toolkit in action:

```bash
cd examples
python basic_example.py
```

This will demonstrate:
1. Single variant annotation
2. Batch processing 
3. CSV file processing
4. Position-based annotation
5. Score interpretation

## Input Formats

### Variant IDs
- **dbSNP rsIDs**: `rs238242`
- **HGVS notation**: `chr1:12345:A>G` 
- **Transcript notation**: `NM_000551.3:c.1521_1523delCTT`

### CSV Format
Your CSV file should have a column with variant IDs:
```csv
variant_id,gene,condition
rs238242,VHL,VHL syndrome
rs7527068,TCF7L2,Type 2 diabetes
chr1:12345:A>G,Unknown,Unknown
```

## Data Sources

| Source | Type | Information Provided |
|--------|------|---------------------|
| MyVariant.info | Variant | CADD scores, ClinVar, gnomAD frequencies |
| Ensembl VEP | Variant | Consequence predictions, gene annotations |
| ClinVar | Variant | Clinical significance, disease associations |
| UCSC PhyloP | Position | Evolutionary conservation scores |
| UCSC PhastCons | Position | Conservation elements |

## Scoring System

The toolkit provides a composite score (0-1) based on:

- **CADD Score**: Higher CADD = higher pathogenicity score
- **ClinVar Significance**: Pathogenic variants get high scores  
- **Population Frequency**: Rarer variants get higher scores
- **Conservation**: More conserved positions get higher scores

### Score Interpretation
- **> 0.6**: HIGH IMPACT - Likely pathogenic
- **0.3-0.6**: MODERATE IMPACT - Possibly significant  
- **0.1-0.3**: LOW IMPACT - Likely benign but notable
- **< 0.1**: MINIMAL IMPACT - Likely benign

## API Reference

### AnnotationPipeline

Main interface for annotation tasks.

#### Methods

- `run_single_variant(variant_id, show_results=True)` - Annotate one variant
- `run_batch_variants(variant_list)` - Process multiple variants
- `run_from_csv(input_file, variant_column, output_file)` - Process CSV file
- `run_position(chromosome, position, ref, alt)` - Annotate genomic position

### GenomicPosition

Represents a genomic position.

```python
from genomic_annotator import GenomicPosition

pos = GenomicPosition("chr1", 12345, "A", "G")
print(pos.to_hgvs())  # chr1:g.12345A>G
```

### Custom Usage

For more control, use the underlying classes:

```python
from genomic_annotator import GenomicAnnotator

annotator = GenomicAnnotator()
results = annotator.annotate_variant("rs238242")
scores = annotator.score_variant(results)
```

## Rate Limiting

All API calls are rate-limited:
- **Default**: 0.6 seconds between requests
- **Automatic**: Built-in retry logic
- **Configurable**: Adjust delays if needed

## Error Handling

- Network timeouts
- API rate limits  
- Invalid variant IDs
- Missing data

Failed annotations are marked as unsuccessful but don't stop batch processing.

## Example Output

```
ANNOTATION RESULTS FOR: rs238242

MYVARIANT:
  Success: True
  cadd: {'phred': 25.3}
  clinvar: {'rcv': [{'clinical_significance': 'Pathogenic'}]}
  gnomad_exome: {'af': 0.0001}

ENSEMBL_VEP:
  Success: True  
  consequence_terms: ['missense_variant']
  gene_symbol: VHL

CLINVAR:
  Success: True
  title: VHL, 5-BP DEL/22-BP INS, NT163

SCORE: 0.847
DETAILS: {'cadd_phred': 25.3, 'clinvar_pathogenic': True, 'gnomad_exome_af': 0.0001}
```

