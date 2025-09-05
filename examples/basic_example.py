#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from genomic_annotator import AnnotationPipeline, GenomicPosition, create_sample_data
import pandas as pd


def example_single_variant():
    """Example 1: Annotate a single variant"""
    print("=" * 60)
    print("EXAMPLE 1: Single Variant Annotation")
    print("=" * 60)
    
    pipeline = AnnotationPipeline()
    
    # Annotate rs238242 (a well-known variant)
    result = pipeline.run_single_variant("rs238242", show_results=True)
    
    print(f"\nFinal score: {result['scores']['total_score']:.3f}")


def example_batch_processing():
    """Example 2: Process multiple variants"""
    print("\n" + "=" * 60)
    print("EXAMPLE 2: Batch Processing")
    print("=" * 60)
    
    pipeline = AnnotationPipeline()
    
    # List of variants to process
    variants = [
        "rs238242",      # VHL gene variant
        "rs7527068",     # TCF7L2 variant  
        "rs142513484",   # BRCA1 variant
        "rs12345"        # Non-existent variant (will fail gracefully)
    ]
    
    print(f"Processing {len(variants)} variants...")
    
    # Process batch
    results_df = pipeline.run_batch_variants(variants)
    
    # Display results
    print("\nResults Summary:")
    print(results_df[["variant_id", "total_score", "successful_annotators", "annotator_names"]].to_string(index=False))
    
    # Get summary stats
    stats = pipeline.get_summary_stats(results_df)
    print(f"\nSummary Statistics:")
    print(f"  Total variants processed: {stats['total_variants']}")
    print(f"  Successfully annotated: {stats['successfully_annotated']}")
    print(f"  Mean score: {stats['mean_score']:.3f}")
    print(f"  High-impact variants (score > 0.5): {stats['high_impact_variants']}")


def example_csv_processing():
    """Example 3: Process variants from CSV file"""
    print("\n" + "=" * 60)
    print("EXAMPLE 3: CSV File Processing") 
    print("=" * 60)
    
    # Create sample data
    variants_df, positions_df = create_sample_data()
    
    # Save sample CSV
    sample_file = "examples/data/sample_variants.csv"
    os.makedirs("examples/data", exist_ok=True)
    variants_df.to_csv(sample_file, index=False)
    print(f"Created sample file: {sample_file}")
    print("Sample data:")
    print(variants_df.to_string(index=False))
    
    # Process CSV
    pipeline = AnnotationPipeline()
    results_df = pipeline.run_from_csv(
        input_file=sample_file,
        variant_column="variant_id",
        output_file="examples/data/annotated_results.csv"
    )
    
    print(f"\nProcessed {len(results_df)} variants")
    print("Results saved to: examples/data/annotated_results.csv")
    
    # Show key columns
    key_cols = ["variant_id", "gene", "total_score", "successful_annotators"]
    available_cols = [col for col in key_cols if col in results_df.columns]
    print("\nKey results:")
    print(results_df[available_cols].to_string(index=False))


def example_position_annotation():
    """Example 4: Annotate genomic positions"""
    print("\n" + "=" * 60) 
    print("EXAMPLE 4: Position-based Annotation")
    print("=" * 60)
    
    pipeline = AnnotationPipeline()
    
    # Annotate a specific genomic position
    result = pipeline.run_position(
        chromosome="chr1",
        position=12345,
        ref="A", 
        alt="G"
    )
    
    print("Position annotation completed!")


def example_custom_scoring():
    """Example 5: Understanding the scoring system"""
    print("\n" + "=" * 60)
    print("EXAMPLE 5: Understanding Scores")
    print("=" * 60)
    
    pipeline = AnnotationPipeline()
    
    # Process a few variants and explain scores
    test_variants = ["rs238242", "rs7527068"]
    
    for variant in test_variants:
        print(f"\n--- Analyzing {variant} ---")
        result = pipeline.run_single_variant(variant, show_results=False)
        
        score = result['scores']['total_score']
        details = result['scores']['details']
        
        print(f"Overall Score: {score:.3f}")
        print("Score breakdown:")
        for key, value in details.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.3f}")
            else:
                print(f"  {key}: {value}")
        
        # Interpretation
        if score > 0.6:
            interpretation = "HIGH IMPACT - Likely pathogenic"
        elif score > 0.3:
            interpretation = "MODERATE IMPACT - Possibly significant"
        elif score > 0.1:
            interpretation = "LOW IMPACT - Likely benign but notable"
        else:
            interpretation = "MINIMAL IMPACT - Likely benign"
        
        print(f"Interpretation: {interpretation}")


if __name__ == "__main__":
    print("Genomic Annotator SDK - Examples")
    print("================================")
    
    try:
        example_single_variant()
        example_batch_processing()
        example_csv_processing()
        example_position_annotation()
        example_custom_scoring()
        
        print("\n" + "=" * 60)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("Make sure you have an internet connection for API calls.")
