# Aim

Sometimes different properties of features need to be examined, for example, the proportion of a genomic region covered by . Overlap features with genomic regions adding an additional column with information depending on the feature class (binary, numeric, integer...).

# Input

Two BED files (genomic regions and features) with at least three cols:

  * Chromosome ID
  * Genomic region start
  * Genomic region end

Some classes demand the input feature BED file to have one aditional column.

# Usage

```{bash}
bash custom_intersect.sh genomic_regions.bed feature.bed class
```

# Options

Specify the feature class. 6 classes are accepted:

* Binary
* Ratio
* Ratio.value
* Integer (count)
* Numeric
* Numeric.sum

# Example

## Input
```{bash}
> genomic_regions.bed
chrA01  1  100
chrA01  101  200
chrA01  201  300
chrA01  301  400
chrA02  1  101

> feature.bed
chrA01  1  5
chrA01  12  14
chrA01  189  192
chrA01  432  433
```

## Execution as ratio

```{bash}
bash custom_intersect.sh genomic_regions.bed feature.bed ratio
```

## Output as ratio

```{bash}
> genomic_regions.bed
chrA01  1  100  0.06
chrA01  101  200  0.03
chrA01  201  300  0.00
chrA01  301  400  0.00
chrA02  1  101  0.00
```

## Execution as integer
```{bash}
bash custom_intersect.sh genomic_regions.bed feature.bed integer
```

## Output as integer

```{bash}
> genomic_regions.bed
chrA01  1  100  6
chrA01  101  200  3
chrA01  201  300  0
chrA01  301  400  0
chrA02  1  101  0
```
