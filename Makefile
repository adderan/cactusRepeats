progressiveCactusDir=${PWD}/../progressiveCactus-repeats/
cactusDir=${progressiveCactusDir}/submodules/cactus
sonLibDir=${progressiveCactusDir}/submodules/sonLib
twopacoDir=${progressiveCactusDir}/submodules/TwoPaCo


#Trimmed sequence
seqStart=3000000
seqStop=3300000
chromosome=chr1

#TwoPaCo
seedLength=19
filterSize=10


lastzXX=${sonLibDir}/bin/cPecanLastz
seedProbabilitiesXX=${cactusDir}/bin/cactus_repeats_seedProbabilities
twopacoXX=${twopacoDir}/build/graphconstructor/twopaco
mafComparatorXX=${PWD}/maftools/bin/mafComparator

sequence=${PWD}/mm10.chr1.fa



benchmarksDir=${PWD}/benchmarks
alignmentsDir=${benchmarksDir}/alignments
seedProbabilitiesDir=${benchmarksDir}/seedProbabilities

seqIndexDir=${benchmarksDir}/seqIndex
graphDir=${benchmarksDir}/graph_output
trimmedSeq=${benchmarksDir}/trimmed_sequence.fa
deBruijnGraph=${benchmarksDir}/de_bruijn.bin
seedCounts=${benchmarksDir}/seed_counts_lastz.txt


#Sampled alignment
seedMultiplicityProbTable=${seedProbabilitiesDir}/seedMultiplicity.txt
seedMultiplicitySamplingStats=${seedProbabilitiesDir}/seedMultiplicity.stats.txt
seedMultiplicityAlignment=${alignmentsDir}/seedMultiplicity.cigar
seedMultiplicityAlignmentStats=${alignmentsDir}/seedMultiplicity.stats.txt

sampledLastzArgs = --format=cigar --notrivial --sampleSeeds ${trimmedSeq}[unmask][multiple] ${trimmedSeq}[unmask][multiple]

maskedLastzArgs = --format=cigar --notrivial ${trimmedSeq}[multiple] ${trimmedSeq}[multiple]

maskedStats=${alignmentsDir}/masked_stats.txt
maskedAlignment=${alignmentsDir}/masked_alignment.cigar

${benchmarksDir}:
	mkdir -p ${benchmarksDir}
${alignmentsDir}:
	mkdir -p ${alignmentsDir}

${seedProbabilitiesDir}:
	mkdir -p ${seedProbabilitiesDir}

${deBruijnGraph}: ${benchmarksDir} ${trimmedSeq} ${twopacoXX}
	mkdir -p ${graphDir}
	${twopacoXX} -f ${filterSize} -k ${seedLength} --tmpdir ${graphDir} -o ${deBruijnGraph} ${trimmedSeq}

${seedCounts}: ${benchmarksDir} ${trimmedSeq} ${lastzXX}
	${lastzXX} --tableonly=count ${trimmedSeq}[unmask][multiple] > ${seedCounts}

${seqIndexDir}: ${sequence}
	rm -rf ${seqIndexDir}
	mkdir -p ${seqIndexDir}
	cd ${seqIndexDir} && samtools faidx ${sequence}

${trimmedSeq}: ${seqIndexDir}
	cd ${seqIndexDir} && samtools faidx ${sequence} ${chromosome}:${seqStart}-${seqStop} > ${trimmedSeq}

${maskedAlignment} ${maskedStats}: ${lastzXX} ${trimmedSeq} ${alignmentsDir}
	${lastzXX} ${maskedLastzArgs} --stats=${maskedStats} > ${maskedAlignment}


#Sampled alignments

${seedMultiplicityProbTable}: ${seedCounts} ${deBruijnGraph} ${trimmedSeq} ${seedProbabilitiesXX} ${seedProbabilitiesDir}
	${seedProbabilitiesXX} --seedMultiplicity --graphFile ${deBruijnGraph} --lastzFile ${seedCounts} --sequenceFile ${trimmedSeq} --statsFile ${seedMultiplicitySamplingStats} > ${seedMultiplicityProbTable}

${seedMultiplicityAlignment} ${seedMultiplicityAlignmentStats}: ${seedMultiplicityProbTable} ${lastzXX} ${trimmedSeq} ${alignmentsDir}
	${lastzXX} ${sampledLastzArgs} --seedSamplingProbabilities=${seedMultiplicityProbTable} --stats=${seedMultiplicityAlignmentStats} > ${seedMultiplicityAlignment}



benchmarks: ${maskedAlignment} ${seedMultiplicityAlignment}

clean:
	rm -rf ${benchmarksDir}


