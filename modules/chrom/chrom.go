package chromosomes

import (
	"math/rand"
)


type Sex int 
const (
	BOY Sex = iota
	GIRL
	NUM_SEXES
)

type Parent int
const (
	DADDY Parent = iota
	MOMMY
	NUM_PARENTS
)

type Arm int
const (
	P Arm = iota
	Q
	NUM_ARMS
)

type ArmData int
const (
	ARMSTART ArmData = iota
	ARMLENGTH
	NUM_ARMDATA
)

var ChromosomeCount int
var GenomeBitCount int // total number of bits in the genome
var GenomArrSize int // number of bytes in the genome array

// Chromosome index, arm: 0=P 1=Q, 0=start 1=length 
var ChromosomeDef [][NUM_ARMS][NUM_ARMDATA]int

func InitChromosomes(count int) {
	ChromosomeCount = count
	ChromosomeDef = make([][NUM_ARMS][NUM_ARMDATA]int, count)
}

func DefineChromosome(idx, pstart, plen, qstart, qlen int) {
	ChromosomeDef[idx][P][ARMSTART] = pstart
	ChromosomeDef[idx][P][ARMLENGTH] = plen
	ChromosomeDef[idx][Q][ARMSTART] = qstart
	ChromosomeDef[idx][Q][ARMLENGTH] = qlen

	// add the combined length of arms to the total length
	// of the genome
	GenomeBitCount += plen + qlen
	GenomeArrSize = *(GenomBitCount + 7) / 8
}

func MakeChildChromosomes(mom, dad [][NUM_PARENTS]uint8, sexOfChild Sex) [][NUM_PARENTS]uint8 {
	genome := make([][NUM_PARENTS]uint8, GenomeArrSize)

	// Gooo goo gah gah
	dadMask := makeMask(sexOfChild, DADDY)
	momMask := makeMask(sexOfChild, MOMMY)
	for arrIdx := 0; arrIdx < GenomeArrSize; ++arrIdx {
		{
			grampaGenes := dad[arrIdx][DADDY] & (^dadMask[arrIdx])
			grammaGenes := dad[arrIdx][MOMMY] & dadMask[arrIdx]
			genome[arrIdx][DADDY] = grampaGenes | grammaGenes
		}
		{
			grampaGenes := mom[arrIdx][DADDY] & (^momMask[arrIdx])
			grammaGenes := mom[arrIdx][MOMMY] & momMask[arrIdx]
			genome[arrIdx][MOMMY] = grampaGenes | grammaGenes
		}
	}

	return genome
}

func makeMask(sexOfChild Sex, parent Parent) []uint8 {
	mask := make([]uint8, GenomeArraySize)
	var maskBitIdx int
	var startChrom int

	// If we're talking about the mommy, the X chromosomes behave like
	// autosomal chromosomes as far as crossovers are concerned, so start
	// the crossover calculations at index 0 which is the sex chromosomes.
	// If it's the daddy however, then we need to do something special with
	// the sex chromosomes so start crossover calculations at index 1 instead.
	if parent == DADDY {
		startChrom = 1

		plen := ChromosomeDef[0][P][ARMLENGTH]
		qlen := ChromosomeDef[0][Q][ARMLENGTH]
		if sexOfChild == GIRL {
			// Eetsagirl, so we need to get X chromosome from daddy,
			// so set all bits to 1's to inherit from daddy's maternal side.
			// Conversely, leaving as 0's would cause the child to inherit
			// daddy's paternal side, so we don't need to flip any bits for
			// boys.
			for bit := maskBitIdx; bit < plen + qlen; ++bit {
				mask[bit / 8] |= 1 << (bit % 8)
			}
		}

		// Advance to the next chromosome
		maskBitIdx = plen + qlen
	}

	// For each regular chromosome, calculate a random crossover index
	// on each arm and a random centromere, then use it to create
	// the mask bits for that chromosome.
	for chrom := startChrom; chrom < ChromosomeCount; ++chrom {
		// grab arm lengths for convenience
		plen := ChromosomeDef[chrom][P][ARMLENGTH]
		qlen := ChromosomeDef[chrom][Q][ARMLENGTH]

		// Calculate genome bit indices of crossovers
		startBitIdx := maskBitIdx
		crossBitIdxP := maskBitIdx + rand.Intn(plen)
		crossBitIdxQ := maskBitIdx + plen + rand.Intn(qlen)
		endBitIdx := maskBitIdx + plen + qlen

		// Choose which centromere
		centromere := rand.Intn(NUM_PARENTS)
		if centromere == DADDY {
			// using daddy's centromere so flip bits on the ends
			for bit := startBitIdx; bit < crossBitIdxP; ++bit {
				mask[bit / 8] |= 1 << (bit % 8)
			}
			for bit := crossBitIdxQ; bit < endBitIdx; ++bit {
				mask[bit / 8] |= 1 << (bit % 8)
			}
		} else {
			// using mommy's centromere so flip bits in the middle
			for bit := crossBitIdxP; bit < crossBitIdxQ; ++bit {
				mask[bit / 8] = 1 << (bit % 8)
			}
		}

		// Advance to the next chromosome
		maskBitIdx += plen + qlen
	}

	return mask
}

