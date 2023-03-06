# Simulation studies

We conducted simulation studies to examine the importance of accounting for effect heterogeneity in fine-mapping and assess the utility of SharePro in GxE analysis.

## Prepare genotype
```
## select individuals
head -50000 ~/scratch/UKB_Geno/1.fam > UKB_A.fam

## sample genotypes
while read a b c d; do sed "s/LOCI/$a/g" LOCI.sh | sed "s/CHR/$b/g" | sed "s/ST/$c/g" | sed "s/ED/$d/g" > $a\.sh; done < random.gtf
for i in $(seq 1 5); do ./Locus$i\.sh; done
```

## Simulate traits

```
./simulate_traits.sh LOCI KC B1 B2
```

## Detect GxE with PLINK

```
./run_plink.sh LOCI KC_B1_B2
```

## Detect GxE with GEM

```
./run_gem.sh LOCI KC_B1_B2
```

## Detect GxE with SharePro

```
./run_SH.sh LOCI KC_B1_B2
```

## Visualize results

```
Rscript ../doc/plot_sharepro_gxe_sim.R
```