"""
Genomic Annotator SDK
"""

__version__ = "1.0.0"

# Main imports
from .annotator import (
    GenomicAnnotator,
    GenomicPosition,
    GenomicRegion,
    AnnotationResult,
    MyVariantAnnotator,
    EnsemblVEPAnnotator,
    ClinVarAnnotator,
    ConservationAnnotator,
    VariantScorer,
    display_results
)

from .pipeline import AnnotationPipeline, create_sample_data

__all__ = [
    "GenomicAnnotator",
    "AnnotationPipeline", 
    "GenomicPosition",
    "GenomicRegion",
    "AnnotationResult",
    "MyVariantAnnotator",
    "EnsemblVEPAnnotator", 
    "ClinVarAnnotator",
    "ConservationAnnotator",
    "VariantScorer",
    "display_results",
    "create_sample_data"
]

