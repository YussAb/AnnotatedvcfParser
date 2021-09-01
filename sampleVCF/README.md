# Test scripts

Some vcf files from different sources to fully test the script

```bash
python3 AnnotatedvcfParser.py -v AnnotatedvcfParser/sampleVCF/annotated_neusomatic_pass.vcf  -o output.csv -t s
```

```bash
python3 AnnotatedvcfParser.py -v nightly-civic_accepted.vcf -o output.csv -t v
```

```bash
python3 AnnotatedvcfParser.py -v Mutect2_filtered_HKNPC-090T_vs_HKNPC-090N_cosmic_filter.vcf -o output.csv -t s,v
```