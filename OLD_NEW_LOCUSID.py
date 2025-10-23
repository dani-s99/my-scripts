#!/usr/bin/env python3
"""
map_locus_tags.py
Mapea locus_tag antiguo -> locus_tag nuevo entre dos archivos GenBank
Salida: correspondencia.tsv con columnas:
old_start, old_end, old_locus, new_start, new_end, new_locus, gene, product, overlap_frac
"""

from Bio import SeqIO
import sys

def extract_features(gbk, feature_types=("CDS", "gene")):
    """
    Devuelve lista de features: cada entrada es dict con keys:
      start (1-based inclusive), end (1-based inclusive), locus, gene, product, length
    """
    feats = []
    for rec in SeqIO.parse(gbk, "genbank"):
        for feat in rec.features:
            if feat.type not in feature_types:
                continue
            qualifiers = feat.qualifiers
            # tomar locus_tag preferentemente; si no existe, saltar
            locus = qualifiers.get("locus_tag", [None])[0]
            # a veces locus_tag está ausente; podrías usar protein_id, gene, etc.
            if locus is None:
                continue
            gene = qualifiers.get("gene", ["NA"])[0]
            product = qualifiers.get("product", ["NA"])[0]
            # Biopython: start es 0-based inclusive, end es 0-based exclusive
            start0 = int(feat.location.start)
            end0 = int(feat.location.end)  # exclusive
            # convertimos a 1-based inclusive
            start = start0 + 1
            end = end0
            length = end - start + 1
            feats.append({
                "start": start,
                "end": end,
                "locus": locus,
                "gene": gene,
                "product": product,
                "length": length
            })
    return feats

def overlap_fraction(a_start, a_end, b_start, b_end):
    """Calcula fracción de solapamiento relativa al menor de los dos tamaños."""
    # solapamiento [max(starts), min(ends)]
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    if e < s:
        return 0.0
    overlap = e - s + 1
    len_a = a_end - a_start + 1
    len_b = b_end - b_start + 1
    # fracción relativa al menor de los dos
    return overlap / min(len_a, len_b)

def map_features(old_feats, new_feats, threshold=0.5):
    """
    Para cada old_feat busca new_feats que solapen por al menos 'threshold' (fracción).
    Devuelve lista de tuplas (old, new, overlap_frac)
    """
    mappings = []
    # indexación simple: iterar todas (O(n^2)). Para genomas bacterianos suele ser OK.
    for o in old_feats:
        matches = []
        for n in new_feats:
            frac = overlap_fraction(o["start"], o["end"], n["start"], n["end"])
            if frac >= threshold:
                matches.append((n, frac))
        # si quieres, elegir la nueva con mayor solapamiento
        if matches:
            # ordenar por fracción desc y añadir top matches
            matches.sort(key=lambda x: x[1], reverse=True)
            for n, frac in matches:
                mappings.append((o, n, frac))
        else:
            # opcional: si no hay matches, puedes reportarlo también con new=None
            mappings.append((o, None, 0.0))
    return mappings

def main(old_gbk, new_gbk, out_tsv="correspondencia.tsv", threshold=0.5):
    old_feats = extract_features(old_gbk)
    new_feats = extract_features(new_gbk)
    mappings = map_features(old_feats, new_feats, threshold=threshold)
    with open(out_tsv, "w") as out:
        out.write("old_start\told_end\told_locus\tnew_start\tnew_end\tnew_locus\tgene\tproduct\toverlap_frac\n")
        for o, n, frac in mappings:
            if n is None:
                out.write(f"{o['start']}\t{o['end']}\t{o['locus']}\tNA\tNA\tNA\t{o['gene']}\t{o['product']}\t{frac:.3f}\n")
            else:
                out.write(f"{o['start']}\t{o['end']}\t{o['locus']}\t{n['start']}\t{n['end']}\t{n['locus']}\t{o['gene']}\t{o['product']}\t{frac:.3f}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: map_locus_tags.py antiguo.gbk nuevo.gbk [out.tsv] [overlap_threshold]")
        sys.exit(1)
    old_gbk = sys.argv[1]
    new_gbk = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else "correspondencia.tsv"
    thr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    main(old_gbk, new_gbk, out_tsv=out, threshold=thr)

