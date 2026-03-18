"""
Variant Annotation via Ensembl VEP REST API
============================================
Annotates GWAS variants with consequence, nearest gene,
rsID, and functional information using the Ensembl VEP API.

Usage:
    from annotation import annotate_variants
    df_annotated = annotate_variants(df_suggestive)
"""

import pandas as pd
import numpy as np
import requests
import time
import json

ENSEMBL_VEP_URL = "https://rest.ensembl.org/vep/human/region"
ENSEMBL_GENE_URL = "https://rest.ensembl.org/overlap/region/human"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
BATCH_SIZE = 200
RATE_LIMIT_SLEEP = 0.1   # seconds between requests


def locus_to_vep_format(row):
    """Convert GWAS locus to Ensembl VEP input format."""
    chrom = row['CHR'].replace('chr', '')
    pos   = int(row['POS'])
    ref   = row['REF']
    alt   = row['ALT']
    # VEP format: chrom pos pos ref/alt
    return f"{chrom} {pos} {pos} {ref}/{alt} 1"


def annotate_variants(df, max_variants=None):
    """
    Annotate GWAS variants using Ensembl VEP REST API.

    Parameters
    ----------
    df : DataFrame with columns CHR, POS, REF, ALT
    max_variants : int, optional limit for testing

    Returns
    -------
    DataFrame with added annotation columns:
      - consequence: most severe consequence
      - gene_name: nearest gene symbol
      - rsid: dbSNP rsID if available
      - biotype: transcript biotype
      - impact: LOW/MODERATE/HIGH/MODIFIER
    """
    if max_variants:
        df = df.head(max_variants).copy()

    print(f"Annotating {len(df):,} variants via Ensembl VEP...")

    results = []
    batches = [df.iloc[i:i+BATCH_SIZE] for i in range(0, len(df), BATCH_SIZE)]

    for batch_idx, batch in enumerate(batches):
        print(f"  Batch {batch_idx+1}/{len(batches)} ({len(batch)} variants)...")

        vep_inputs = [locus_to_vep_format(row) for _, row in batch.iterrows()]
        payload = {"variants": vep_inputs}

        try:
            response = requests.post(
                ENSEMBL_VEP_URL,
                headers=HEADERS,
                data=json.dumps(payload),
                timeout=60
            )
            response.raise_for_status()
            vep_results = response.json()

        except requests.exceptions.RequestException as e:
            print(f"    API error on batch {batch_idx}: {e}")
            for _, row in batch.iterrows():
                results.append(_empty_annotation(row))
            continue

        # Parse VEP results
        vep_by_input = {r.get('input', ''): r for r in vep_results}

        for _, row in batch.iterrows():
            input_str = locus_to_vep_format(row)
            vep_hit   = vep_by_input.get(input_str, {})
            results.append(_parse_vep_result(row, vep_hit))

        time.sleep(RATE_LIMIT_SLEEP)

    annotations = pd.DataFrame(results)
    df_annotated = df.merge(
        annotations[['locus', 'consequence', 'gene_name', 'rsid',
                     'biotype', 'impact']],
        on='locus', how='left'
    )
    print(f"  ✓ Annotation complete")
    return df_annotated


def _parse_vep_result(row, vep_hit):
    """Extract key fields from a VEP result entry."""
    ann = {
        'locus':       row.get('locus', f"chr{row['CHR']}:{row['POS']}"),
        'consequence': 'intergenic_variant',
        'gene_name':   None,
        'rsid':        None,
        'biotype':     None,
        'impact':      'MODIFIER'
    }

    if not vep_hit:
        return ann

    # rsID
    ann['rsid'] = vep_hit.get('id')

    # Most severe consequence
    ann['consequence'] = vep_hit.get('most_severe_consequence', 'intergenic_variant')

    # Transcript consequences
    tc = vep_hit.get('transcript_consequences', [])
    if tc:
        # Prefer protein-coding transcript
        protein_coding = [t for t in tc if t.get('biotype') == 'protein_coding']
        best = protein_coding[0] if protein_coding else tc[0]
        ann['gene_name'] = best.get('gene_symbol')
        ann['biotype']   = best.get('biotype')
        ann['impact']    = best.get('impact', 'MODIFIER')
    else:
        # Fall back to intergenic annotation
        ig = vep_hit.get('intergenic_consequences', [{}])
        if ig:
            ann['impact'] = ig[0].get('impact', 'MODIFIER')

    return ann


def _empty_annotation(row):
    """Return empty annotation for failed API calls."""
    return {
        'locus':       row.get('locus', f"chr{row['CHR']}:{row['POS']}"),
        'consequence': None,
        'gene_name':   None,
        'rsid':        None,
        'biotype':     None,
        'impact':      None
    }


def get_genes_in_region(chrom, start, end, biotype='protein_coding'):
    """
    Query Ensembl for genes in a genomic region.

    Parameters
    ----------
    chrom : str, e.g. '15' or 'chr15'
    start : int
    end   : int
    biotype : str, default 'protein_coding'

    Returns
    -------
    List of dicts with gene name, coordinates, distance to center
    """
    chrom = chrom.replace('chr', '')
    center = (start + end) // 2

    url = f"{ENSEMBL_GENE_URL}/{chrom}:{start}-{end}"
    params = {"content-type": "application/json", "feature": "gene"}

    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        genes = r.json()
    except Exception as e:
        print(f"Gene query error: {e}")
        return []

    results = []
    for g in genes:
        if biotype and g.get('biotype') != biotype:
            continue
        dist = min(abs(g['start'] - center), abs(g['end'] - center))
        results.append({
            'gene_name': g.get('external_name', '?'),
            'biotype':   g.get('biotype', ''),
            'start':     g['start'],
            'end':       g['end'],
            'strand':    g['strand'],
            'distance_to_center': dist
        })

    return sorted(results, key=lambda x: x['distance_to_center'])
