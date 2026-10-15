# Toxo Genotyper

*Sanger MLST allele calls, electronic PCR-RFLP and ToxoDB genotypes*

A free, browser-based tool for genotyping *Toxoplasma gondii* from Sanger reads of the standard typing markers. It assembles forward and reverse reads and calls the allele at each marker, replacing the manual Geneious Prime workflow. It also cuts each sequence in silico to give the PCR-RFLP type, and matches the multilocus profile to a ToxoDB genotype number. No software installation or computational skills are needed.

**Status:** proof of concept. PCR primers are built in for 13 markers. Reference alleles are built in for all 13: B1 (from published B1 sequences) and 11 others taken from 15 *T. gondii* genomes in ToxoDB (17 for c22-8, which adds CtCo5 and TgCATBr5; 14 for Apico and 13 for CS3, see below) (ARI, BR9, CAST, COUG, DOM2, FOU, GT1, MAS, ME49, P89, PRC2, RH88, RUB, VAND, VEG). PK1 also includes the lab's own Type I, II, III, X and COUG references.

## Reference panel

| Marker | Expected product | Compared region | Distinct alleles among the 15 genomes |
|---|---|---|---|
| PK1 | 902–903 bp | 845 bp | 10 |
| SAG1 | 476 bp | 435 bp | 6 |
| alt.SAG2 | 726–729 bp | 686 bp | 7 |
| 3′-SAG2 | 326–327 bp | 287 bp | 5 |
| SAG3 | 311 bp | 271 bp | 5 |
| BTUB | 411 bp | 371 bp | 5 |
| GRA6 | 344 bp | 306 bp | 9 |
| c22-8 | 521 bp | 485 bp | 11 (of 17 genomes) |
| B1 (multi-copy) | — (inner product ≈ 470 bp) | 439 bp | 12 named alleles from Zhu et al. 2023 |
| c29-2 | 446 bp | 410 bp | 9 |
| L358 | 417 bp (published 418) | 381 bp | 6 |
| Apico | 639 bp (published 640) | 597 bp | 4 (14 strains; ARI, COUG and DOM2 excluded because their available sequence does not span both primers) |
| CS3 | 794 bp outer, 557 bp inner (Pena et al. 2008) | 754 bp | 9 (13 strains; ARI and FOU excluded, assembly gaps over the amplicon) |

Every primer pair was checked against all 15 genomes (both sites found, at most one mismatch). Notes:

- **SAG1:** the forward primer sits in a 27 bp tandem repeat upstream of the start codon and can bind more than one copy. The tool uses the innermost pair of sites, which gives the same 476 bp product in every strain.
- **alt.SAG2:** the built-in pair is the outer PCR. The inner (nested) pair from Su et al. 2006 gives 546 bp and is listed as a note. If the nested product is what gets sequenced, switch to the inner pair.
- **c22-8:** both primers were confirmed in a flanked ME49 sequence (521 bp product). The other strains' references are 491–492 bp genome regions that span the stretch between the primers plus a few bases of each primer site.
- **B1** is a multi-copy gene (about 35 copies), so double peaks are part of the allele: Types II/III carry C/T at position 366 and C/G at 504. For B1 the tool matches double peaks exactly instead of reading them as mixed infection, and tolerates a reference double peak that a sample shows only partly (a minor copy not detected). References are the named alleles of Zhu et al. 2023 (PLoS NTD 17:e0011829, Table 3): GenBank KM243024, MK988572, KM243022 and OQ850747–OQ850751, plus Type I, Type II/III, Type II/III variant A and Atypical 2 reconstructed from the table's positions. The reconstructed ones should be replaced with the lab's own sequences. `test/b1_table3.js` checks the calls against Table 3 (11/11 rows with definite bases agree).
- **L358:** some genome assemblies have gaps (runs of N) near this locus; the tool never lets N count as a primer match. References are trimmed to the amplicon plus 500 bp on each side.
- **PK1 and BTUB** references are ToxoDB transcript or gene sequences without flanks; both primer sites lie inside them, and the products match the published sizes.

**Nested PCR:** Settings has a "Reads come from the inner (nested) PCR" option. For SAG1 and alt.SAG2, which have known inner primers, it compares between the inner primers instead of the outer ones.

**Multi-sequence FASTA files** (for example a GenBank download) are read record by record: the sample name comes from the record title (the "isolate …" part when present), and records without a forward/reverse label are treated as finished sequences, one row each.

**Test against published data:** `test/study_battisti2018.js` re-analyses the GenBank sequences of Battisti et al. 2018 (Vet Parasitol 253:43–47; MG587956–MG588017) and compares the calls with the paper's Table 2. Download the sequences to `data/study/seqs.fasta` first.

When several references are identical over the compared region, the call shows the first name in panel order, and the others are listed as identical ("Type I = 5 more").

## How to use it

1. Open the tool's web page (or download `index.html` and double-click it).
2. Drop in your `.ab1` files, a folder, or a `.zip`. Forward and reverse reads are paired by sample name.
3. Clean samples are called automatically. Samples marked amber or red come with an explanation and the chromatogram at each flagged position.
4. Below the allele table, the electronic PCR-RFLP table gives the RFLP genotype; select a call to see its virtual gel.
5. Download the allele table (CSV, opens in Excel), the RFLP table, the details and the consensus sequences.

Your sequence files are read inside your own browser and are never uploaded anywhere.

## What it does

| Step | Method |
|---|---|
| Match each read to a marker and orient it | Shared 12-mers with the reference panel |
| Quality trimming (.ab1 only) | Mott's algorithm, error limit 0.05 (the Geneious default) |
| Mixed / heterozygous positions (.ab1 only) | Second peak ≥ 0.33 × main peak, written as an IUPAC code |
| Forward + reverse assembly | Overlap alignment with affine gaps; disagreements are excluded and reported |
| Comparison with references | Alignment to every reference over the region between the marker's PCR primers (or, without primers or with the full-length setting, the region all references cover). References can be amplicons or whole genes |
| Allele call | Zero differences required. Anything else is flagged for review, with an NCBI BLAST link |

## Electronic PCR-RFLP

The tool also gives a PCR-RFLP genotype, so that sequence results can be compared with the many studies that report only RFLP types. Each consensus is cut in silico with the enzymes of the multilocus PCR-RFLP scheme on the **nested** PCR product (Su et al. 2006, *Int J Parasitol* 36:841; Su et al. 2010, *Parasitology* 137:1, Table 2). A virtual gel is shown for each sample, next to the Type I, II and III patterns.

| Locus | Enzymes | Nested product | Typed from marker |
|---|---|---|---|
| SAG1 | Sau96I + HaeII | 390 bp | SAG1 |
| 5′-SAG2 | MboI | 242 bp | alt.SAG2 (inside the outer product) |
| 3′-SAG2 | HhaI | 222 bp | 3′-SAG2 |
| alt.SAG2 | HinfI + TaqI | 546 bp | alt.SAG2 |
| SAG3 | NciI | 225 bp | SAG3 |
| BTUB | BsiEI + TaqI | 411 bp | BTUB |
| GRA6 | MseI | 344 bp | GRA6 |
| c22-8 | BsmAI + MboII | 521 bp | c22-8 |
| c29-2 | HpyCH4IV + RsaI | 446 bp | c29-2 |
| L358 | HaeIII + NlaIII | 418 bp | L358 |
| PK1 | AvaI + RsaI | 903 bp | PK1 |
| Apico | AflII + DdeI | 640 bp | Apico |
| CS3 | NlaIII + MboI | 557 bp | CS3 (Pena et al. 2008) |

How a pattern is named:

- The panel references are cut the same way and grouped by banding pattern. Two patterns count as the same when they have the same visible bands (≥ 30 bp), each within about 5%, which is what a 2.5–3% agarose gel can separate.
- A group is named after the archetype strains in it (I = GT1, II = ME49, III = VEG). "I or II" means the enzymes at that locus do not separate those types; SAG1 cannot separate II and III, as the paper notes.
- Non-archetypal groups take the published name (u-1, u-2) when a panel strain carries it in Su et al. Otherwise they are "new".
- SAG2 (5′+3′) combines the two halves as in Su et al. 2010 Table 3.
- Where the consensus does not reach the end of the nested product, or has unread bases, those bases are taken from the closest reference. If the pattern there differs between references, the locus is reported as ND (not determined).
- A cut site that depends on a double peak gives two possible patterns. The site is flagged, because this can mean a mixed infection.

**ToxoDB genotype.** The multilocus e-RFLP profile is matched against the ToxoDB PCR-RFLP genotype reference table (190 genotypes, stored in `src/rflp.json`). Only 5 genotypes in the table have a CS3 type (#8, #11, #13, #108, #232); for the rest CS3 is not used. The tool reports every genotype whose profile fits all loci typed. With all 11 loci typed this is usually one genotype; with only a few loci typed, several can fit and the tool lists them. In the table's SAG1 column, II stands for "II or III". A pattern not seen in the panel ("new") is allowed to match any u- allele. When nothing fits, the closest genotypes and the loci where they differ are given. Typing the panel's own genome strains this way gives GT1 and RH88 #10, ME49 #1, VEG #2, MAS #17, COUG #66, FOU #6, TgCatPRC2 #9, P89 #8 (its CS3 pattern, III, also matches the CS3 type listed for #8), VAND #60, BR9 #42, CAST #28, RUB #98, ARI #5 and DOM2 #203 (ARI, COUG and DOM2 have no Apico reference).

Validation (`node test/rflp_published.js`): the panel's genome strains were typed and compared with the published PCR-RFLP types of the same strains (GT1, RH88, VEG, CAST, COUG, MAS, VAND, TgCatBr5, TgCatBr9; Su et al. 2006 Table 2, Su et al. 2010 Table 3). 79 of 81 agree. Both disagreements are CAST (SAG3 and Apico), where the ToxoDB CAST genome sequence itself carries a different allele from the CAST isolate typed in 2006. `node test/rflp_samples.js <folder or fasta> [nested]` prints the e-RFLP of real samples.

## Reference panel format

Load your own references from the **Reference panel** section of the page. The file is FASTA, with headers in this form:

```
>MARKER|Allele name|optional note
>SAG3|Type II|GenBank AB123456
```

Primers go on a header line with no sequence below it:

```
>SAG3|primers|FORWARD_PRIMER|REVERSE_PRIMER
```

## Files

- `index.html`: the complete tool in one file
- `src/core.js`: analysis code (also used by the tests)
- `src/app.html`: page template; `build.py` combines it with `core.js`, the panel and the primers to produce `index.html`
- `src/panel.json`: built-in reference alleles
- `src/primers.json`: built-in PCR primers and their sources
- `src/rflp.js`, `src/rflp.json`: electronic PCR-RFLP (enzymes, nested primers, published types used for naming and testing)
- `add_marker.py`, `check_marker.js`: add a marker's ToxoDB sequences to the panel, then check primer sites, alleles and simulated calls (`python3 add_marker.py GRA6 gra6.fasta 'genomic sequence'`, then `node check_marker.js GRA6` and `python3 build.py`)
- `test/`: Node test runner and a synthetic `.ab1` generator

## Credits

Developed by Omar S. Harb (VEuPathDB, University of Pennsylvania) for the *Toxoplasma* research community. The code was written with assistance from Claude, an AI model made by Anthropic; the design, reference data and validation were directed and checked by the author.

PCR primers, restriction enzymes and published genotypes come from the studies cited above and in `src/primers.json` and `src/rflp.json`. Reference sequences are from ToxoDB and GenBank. The ToxoDB genotype numbers come from the ToxoDB PCR-RFLP genotype reference table.
