"""
Genomic Annotator - Multi-level variant annotation system
"""

import requests
import time
import logging
from typing import Dict, Optional, List, Any, Union
from dataclasses import dataclass
from enum import Enum

# Configuration
RATE_LIMIT_DELAY = 0.6
TIMEOUT = 30

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("genomic_annotator")

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "GenomicAnnotator/1.0"})


class AnnotationType(Enum):
    VARIANT = "variant"
    POSITION = "position" 
    REGION = "region"


@dataclass
class GenomicPosition:
    chromosome: str
    position: int
    reference: str = ""
    alternate: str = ""

    def to_hgvs(self) -> str:
        if self.reference and self.alternate:
            return f"{self.chromosome}:g.{self.position}{self.reference}>{self.alternate}"
        return f"{self.chromosome}:g.{self.position}"


@dataclass
class GenomicRegion:
    chromosome: str
    start: int
    end: int

    def to_ucsc_format(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}"


@dataclass
class AnnotationResult:
    source: str
    data: Dict[str, Any]
    success: bool
    error: Optional[str] = None


class BaseAnnotator:
    """Base class for all annotators"""
    
    def __init__(self, name: str, annotation_type: AnnotationType):
        self.name = name
        self.annotation_type = annotation_type
        self._last_request = 0

    def _rate_limit(self):
        elapsed = time.time() - self._last_request
        if elapsed < RATE_LIMIT_DELAY:
            time.sleep(RATE_LIMIT_DELAY - elapsed)
        self._last_request = time.time()

    def _request(self, url: str, params: Dict = None) -> Optional[Dict]:
        try:
            self._rate_limit()
            resp = SESSION.get(url, params=params, timeout=TIMEOUT)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            log.warning(f"Request failed for {url}: {e}")
            return None


class MyVariantAnnotator(BaseAnnotator):
    """MyVariant.info annotator"""
    
    def __init__(self):
        super().__init__("myvariant", AnnotationType.VARIANT)
        self.base_url = "https://myvariant.info/v1/variant/"
        self.fields = "cadd.phred,clinvar,dbnsfp,gnomad_exome.af,gnomad_genome.af"

    def annotate(self, variant_id: str) -> AnnotationResult:
        url = f"{self.base_url}{variant_id}"
        params = {"fields": self.fields}
        data = self._request(url, params)
        
        success = data is not None and not (isinstance(data, dict) and 
                                          all(k.startswith('_') for k in data.keys()))
        
        return AnnotationResult(
            source=self.name,
            data=data or {},
            success=success,
            error=None if success else "No data found"
        )


class EnsemblVEPAnnotator(BaseAnnotator):
    """Ensembl VEP annotator"""
    
    def __init__(self):
        super().__init__("ensembl_vep", AnnotationType.VARIANT)
        self.base_url = "https://rest.ensembl.org/vep/human/id/"

    def annotate(self, variant_id: str) -> AnnotationResult:
        url = f"{self.base_url}{variant_id}"
        params = {"content-type": "application/json"}
        data = self._request(url, params)
        
        success = data is not None and len(data) > 0 if isinstance(data, list) else data is not None
        if success and isinstance(data, list):
            data = data[0]

        return AnnotationResult(
            source=self.name,
            data=data or {},
            success=success,
            error=None if success else "No data found"
        )


class ClinVarAnnotator(BaseAnnotator):
    """ClinVar annotator"""
    
    def __init__(self):
        super().__init__("clinvar", AnnotationType.VARIANT)
        self.eutils_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def annotate(self, variant_id: str) -> AnnotationResult:
        # Search for variant
        search_url = f"{self.eutils_base}esearch.fcgi"
        search_params = {"db": "clinvar", "term": variant_id, "retmode": "json", "retmax": 5}
        search_data = self._request(search_url, search_params)
        
        if not search_data or not search_data.get("esearchresult", {}).get("idlist"):
            return AnnotationResult(
                source=self.name, data={}, success=False, error="Variant not found in ClinVar"
            )

        # Get summary
        clinvar_id = search_data["esearchresult"]["idlist"][0]
        summary_url = f"{self.eutils_base}esummary.fcgi"
        summary_params = {"db": "clinvar", "id": clinvar_id, "retmode": "json"}
        summary_data = self._request(summary_url, summary_params)
        
        data = summary_data.get("result", {}).get(clinvar_id, {}) if summary_data else {}
        
        return AnnotationResult(
            source=self.name,
            data=data,
            success=bool(data),
            error=None if data else "No ClinVar data found"
        )


class ConservationAnnotator(BaseAnnotator):
    """UCSC conservation scores"""
    
    def __init__(self, track_name: str = "phyloP100way"):
        super().__init__(f"ucsc_{track_name}", AnnotationType.POSITION)
        self.base_url = "https://api.genome.ucsc.edu/getData/track"
        self.track = track_name

    def annotate(self, position: GenomicPosition) -> AnnotationResult:
        params = {
            "genome": "hg38",
            "track": self.track,
            "chrom": position.chromosome,
            "start": position.position - 1,
            "end": position.position
        }
        data = self._request(self.base_url, params)
        
        return AnnotationResult(
            source=self.name,
            data=data or {},
            success=bool(data),
            error=None if data else f"No {self.track} data found"
        )


class VariantScorer:
    """Simple scoring system for variants"""
    
    def score_variant(self, results: Dict[str, AnnotationResult]) -> Dict[str, Any]:
        total_score = 0.0
        details = {}
        
        # MyVariant scoring
        if "myvariant" in results and results["myvariant"].success:
            data = results["myvariant"].data
            score = 0.0
            
            # CADD score
            if "cadd" in data and isinstance(data["cadd"], dict):
                cadd = self._safe_float(data["cadd"].get("phred"))
                if cadd > 0:
                    details["cadd_phred"] = cadd
                    score += min(cadd / 30.0, 0.4)  # Cap at 0.4
            
            # ClinVar significance
            if "clinvar" in data and isinstance(data["clinvar"], dict):
                rcv = data["clinvar"].get("rcv", [])
                if isinstance(rcv, list):
                    for r in rcv:
                        if isinstance(r, dict):
                            sig = str(r.get("clinical_significance", "")).lower()
                            if "pathogenic" in sig:
                                score += 0.4
                                details["clinvar_pathogenic"] = True
                                break
            
            # Frequency (rarer = higher score)
            for db in ["gnomad_exome", "gnomad_genome"]:
                if db in data and isinstance(data[db], dict):
                    af = self._safe_float(data[db].get("af"))
                    if af is not None:
                        details[f"{db}_af"] = af
                        if af == 0:
                            score += 0.2
                        elif af < 0.001:
                            score += 0.15
                        elif af < 0.01:
                            score += 0.1
                        break
            
            total_score += score
            details["myvariant_score"] = score
        
        return {"total_score": min(total_score, 1.0), "details": details}
    
    def _safe_float(self, value, default=None):
        try:
            if value is None:
                return default
            if isinstance(value, (list, tuple)) and value:
                return float(value[0])
            return float(value)
        except (ValueError, TypeError, IndexError):
            return default


class GenomicAnnotator:
    """Main annotation pipeline"""
    
    def __init__(self):
        self.variant_annotators = {
            "myvariant": MyVariantAnnotator(),
            "ensembl_vep": EnsemblVEPAnnotator(),
            "clinvar": ClinVarAnnotator(),
        }
        self.position_annotators = {
            "phylop": ConservationAnnotator("phyloP100way"),
            "phastcons": ConservationAnnotator("phastCons100way"),
        }
        self.scorer = VariantScorer()
    
    def annotate_variant(self, variant_id: str) -> Dict[str, AnnotationResult]:
        """Annotate a single variant"""
        results = {}
        for name, annotator in self.variant_annotators.items():
            log.info(f"Annotating {variant_id} with {name}")
            results[name] = annotator.annotate(variant_id)
        return results
    
    def annotate_position(self, position: GenomicPosition) -> Dict[str, AnnotationResult]:
        """Annotate a genomic position"""
        results = {}
        for name, annotator in self.position_annotators.items():
            log.info(f"Annotating {position.to_hgvs()} with {name}")
            results[name] = annotator.annotate(position)
        return results
    
    def score_variant(self, results: Dict[str, AnnotationResult]) -> Dict[str, Any]:
        """Score variant results"""
        return self.scorer.score_variant(results)
    
    def annotate_variants_batch(self, variant_ids: List[str]) -> List[Dict[str, Any]]:
        """Batch annotate multiple variants"""
        all_results = []
        
        for variant_id in variant_ids:
            log.info(f"Processing variant: {variant_id}")
            
            # Annotate
            results = self.annotate_variant(variant_id)
            
            # Score
            scores = self.score_variant(results)
            
            # Compile result
            variant_result = {
                "variant_id": variant_id,
                "successful_annotators": [name for name, r in results.items() if r.success],
                "total_score": scores["total_score"],
                "score_details": scores["details"],
                "annotations": {name: r.data for name, r in results.items() if r.success}
            }
            
            all_results.append(variant_result)
        
        return all_results


def display_results(results: Dict[str, AnnotationResult], max_items: int = 3):
    """Display annotation results in a readable format"""
    for source, result in results.items():
        print(f"\n{source.upper()}:")
        print(f"  Success: {result.success}")
        
        if result.success and result.data:
            # Show first few items of data
            items = list(result.data.items())[:max_items]
            for key, value in items:
                if isinstance(value, dict):
                    print(f"  {key}: {dict(list(value.items())[:2])}...")
                elif isinstance(value, list):
                    print(f"  {key}: [{len(value)} items]")
                else:
                    print(f"  {key}: {value}")
            if len(result.data) > max_items:
                print(f"  ... and {len(result.data) - max_items} more fields")
        
        if result.error:
            print(f"  Error: {result.error}")
