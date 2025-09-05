
import pandas as pd
import logging
from typing import List, Dict, Any, Optional
from .annotator import GenomicAnnotator, GenomicPosition, display_results

log = logging.getLogger(__name__)


class AnnotationPipeline:
    
    def __init__(self):
        self.annotator = GenomicAnnotator()
    
    def run_single_variant(self, variant_id: str, show_results: bool = True) -> Dict[str, Any]:
        
        results = self.annotator.annotate_variant(variant_id)
        scores = self.annotator.score_variant(results)
        
        if show_results:
            print(f"\n=== ANNOTATION RESULTS FOR: {variant_id} ===")
            display_results(results)
            print(f"\nSCORE: {scores['total_score']:.3f}")
            if scores['details']:
                print(f"DETAILS: {scores['details']}")
        
        return {
            "variant_id": variant_id,
            "results": results,
            "scores": scores
        }
    
    def run_batch_variants(self, variant_ids: List[str]) -> pd.DataFrame:
        
        log.info(f"Processing {len(variant_ids)} variants...")
        results = self.annotator.annotate_variants_batch(variant_ids)
        
        # Convert to DataFrame
        df_data = []
        for result in results:
            row = {
                "variant_id": result["variant_id"],
                "total_score": result["total_score"],
                "successful_annotators": len(result["successful_annotators"]),
                "annotator_names": ",".join(result["successful_annotators"]),
            }
            
            for key, value in result["score_details"].items():
                row[f"score_{key}"] = value
            
          
            annotations = result["annotations"]
            if "myvariant" in annotations:
                mv_data = annotations["myvariant"]
                if "cadd" in mv_data and isinstance(mv_data["cadd"], dict):
                    row["cadd_phred"] = mv_data["cadd"].get("phred")
                
                if "gnomad_exome" in mv_data and isinstance(mv_data["gnomad_exome"], dict):
                    row["gnomad_af"] = mv_data["gnomad_exome"].get("af")
            
            df_data.append(row)
        
        return pd.DataFrame(df_data)
    
    def run_from_csv(
        self, 
        input_file: str, 
        variant_column: str = "variant_id",
        output_file: Optional[str] = None
    ) -> pd.DataFrame:
        """Process variants from CSV file"""
        
        
        log.info(f"Loading variants from {input_file}")
        df_input = pd.read_csv(input_file)
        
        if variant_column not in df_input.columns:
            raise ValueError(f"Column '{variant_column}' not found in {input_file}")
        
       
        variants = df_input[variant_column].tolist()
        
       
        df_results = self.run_batch_variants(variants)
        
     
        df_final = df_input.merge(df_results, on="variant_id", how="left")
        
        
        if output_file:
            log.info(f"Saving results to {output_file}")
            df_final.to_csv(output_file, index=False)
        
        return df_final
    
    def run_position(self, chromosome: str, position: int, ref: str = "", alt: str = "") -> Dict[str, Any]:
     
        
        pos = GenomicPosition(chromosome=chromosome, position=position, reference=ref, alternate=alt)
        results = self.annotator.annotate_position(pos)
        
        print(f"\n=== POSITION ANNOTATION: {pos.to_hgvs()} ===")
        display_results(results)
        
        return {"position": pos, "results": results}
    
    def get_summary_stats(self, df: pd.DataFrame) -> Dict[str, Any]:
      
        
        stats = {
            "total_variants": len(df),
            "successfully_annotated": len(df[df["successful_annotators"] > 0]),
            "mean_score": df["total_score"].mean(),
            "high_impact_variants": len(df[df["total_score"] > 0.5]),
            "annotator_success_rates": {}
        }
        
        # Calculate per-annotator success rates
        for annotator in ["myvariant", "ensembl_vep", "clinvar"]:
            col_name = f"{annotator}_success"
            if col_name in df.columns:
                stats["annotator_success_rates"][annotator] = df[col_name].sum() / len(df)
        
        return stats


def create_sample_data():
    
    # Sample variants
    variants_df = pd.DataFrame({
        "variant_id": [
            "rs238242",
            "rs7527068", 
            "rs142513484",
            "chr1:12345:A>G",
            "chr2:67890:C>T"
        ],
        "gene": ["VHL", "TCF7L2", "BRCA1", "Unknown", "Unknown"],
        "condition": ["VHL syndrome", "Diabetes", "Breast cancer", "Unknown", "Unknown"]
    })
    
    # Sample positions
    positions_df = pd.DataFrame({
        "chromosome": ["chr1", "chr2", "chr3", "chrX"],
        "position": [12345, 67890, 111222, 54321],
        "reference": ["A", "C", "G", "T"],
        "alternate": ["G", "T", "A", "C"]
    })
    
    return variants_df, positions_df


if __name__ == "__main__":
    # Example usage
    pipeline = AnnotationPipeline()
    
    # Test single variant
    pipeline.run_single_variant("rs238242")
    
    # Test batch
    variants = ["rs238242", "rs7527068"]
    results_df = pipeline.run_batch_variants(variants)
    print("\nBatch Results:")
    print(results_df)
