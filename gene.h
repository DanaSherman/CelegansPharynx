#ifndef GENE_H_
#define GENE_H_

#include "base.h"
#include "params.h"

struct GeneRange {
	double min;
	double step;
	double max;
	double known;
};

extern const std::array<GeneRange, Params::NUM_GENES> GENE_RANGES;

#endif	// GENE_H_
