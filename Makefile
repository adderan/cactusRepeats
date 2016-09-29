progressiveCactusDir=${PWD}/../progressiveCactus-repeats/
cactusDir=${progressiveCactusDir}/submodules/cactus
sonLibDir=${progressiveCactusDir}/submodules/sonLib
twopacoDir=${progressiveCactusDir}/submodules/TwoPaCo


#Trimmed sequence
seqStart=3000000
seqStop=3300000
chromosome=chr1

globalSeqStart=3000000
globalSeqStop=4000000

#TwoPaCo
seedLength=19
filterSize=10


lastzXX=${sonLibDir}/bin/cPecanLastz
seedProbabilitiesXX=${cactusDir}/bin/cactus_repeats_seedProbabilities
twopacoXX=${twopacoDir}/build/graphconstructor/twopaco
mafComparatorXX=${PWD}/maftools/bin/mafComparator

sequence=${PWD}/mm10.fa



benchmarksDir=${PWD}/benchmarks
alignmentsDir=${benchmarksDir}/alignments
seedProbabilitiesDir=${benchmarksDir}/seedProbabilities

seqIndexDir=${benchmarksDir}/seqIndex
graphDir=${benchmarksDir}/graph_output
trimmedSeq=${benchmarksDir}/trimmed_sequence.fa
globalSeq=${benchmarksDir}/global_sequence.fa
deBruijnGraph=${benchmarksDir}/de_bruijn.bin
seedCounts=${benchmarksDir}/seed_counts_lastz.txt


#Sampled alignment
seedMultiplicityProbTable=${seedProbabilitiesDir}/seedMultiplicity.txt
seedMultiplicitySamplingStats=${seedProbabilitiesDir}/seedMultiplicity.stats.txt
seedMultiplicityAlignment=${alignmentsDir}/seedMultiplicity.maf
seedMultiplicityAlignmentStats=${alignmentsDir}/seedMultiplicity.stats.txt
seedMultiplicityQuality=${alignmentsDir}/seedMultiplicityPrecisionRecall.xml

pathMultiplicityProbTable=${seedProbabilitiesDir}/pathMultiplicity.txt
pathMultiplicitySamplingStats=${seedProbabilitiesDir}/pathMultiplicity.stats.txt
pathMultiplicityAlignment=${alignmentsDir}/pathMultiplicity.maf
pathMultiplicityAlignmentStats=${alignmentsDir}/pathMultiplicity.stats.txt
pathMultiplicityQuality=${alignmentsDir}/pathMultiplicityPrecisionRecall.xml

sampledLastzArgs = --format=maf --notrivial ${trimmedSeq}[unmask][multiple] ${trimmedSeq}[unmask][multiple]

maskedLastzArgs = --format=maf --notrivial ${trimmedSeq}[multiple] ${trimmedSeq}[multiple]
unmaskedLastzArgs = --format=maf --notrivial ${trimmedSeq}[unmask][multiple] ${trimmedSeq}[unmask][multiple]

maskedStats=${alignmentsDir}/masked_stats.txt
maskedAlignment=${alignmentsDir}/masked_alignment.maf

unmaskedStats=${alignmentsDir}/unmasked_stats.txt
unmaskedAlignment=${alignmentsDir}/unmaskedAlignment.maf
unmaskedQuality=${alignmentsDir}/unmaskedPrecisionRecall.xml

${benchmarksDir}:
	mkdir -p ${benchmarksDir}
${alignmentsDir}:
	mkdir -p ${alignmentsDir}

${seedProbabilitiesDir}:
	mkdir -p ${seedProbabilitiesDir}

${deBruijnGraph}: ${benchmarksDir} ${globalSeq} ${twopacoXX}
	mkdir -p ${graphDir}
	${twopacoXX} -f ${filterSize} -k ${seedLength} --tmpdir ${graphDir} -o ${deBruijnGraph} ${globalSeq}

${seedCounts}: ${benchmarksDir} ${globalSeq} ${lastzXX}
	${lastzXX} --tableonly=count ${globalSeq}[unmask][multiple] > ${seedCounts}

${seqIndexDir}: ${sequence}
	rm -rf ${seqIndexDir}
	mkdir -p ${seqIndexDir}
	cd ${seqIndexDir} && samtools faidx ${sequence}

${trimmedSeq}: ${seqIndexDir}
	cd ${seqIndexDir} && samtools faidx ${sequence} ${chromosome}:${seqStart}-${seqStop} > ${trimmedSeq}

${globalSeq}: ${seqIndexDir}
	cd ${seqIndexDir} && samtools faidx ${sequence} ${chromosome}:${globalSeqStart}-${globalSeqStop} > ${globalSeq}

${maskedAlignment} ${maskedStats}: ${lastzXX} ${trimmedSeq} ${alignmentsDir}
	${lastzXX} ${maskedLastzArgs} --stats=${maskedStats} > ${maskedAlignment}

${unmaskedAlignment} ${unmaskedStats} : ${lastzXX} ${trimmedSeq} ${alignmentsDir}
	${lastzXX} ${unmaskedLastzArgs} --stats=${unmaskedStats} > ${unmaskedAlignment}

${unmaskedQuality} : ${unmaskedAlignment} ${maskedAlignment}
	${mafComparatorXX} --maf1=${maskedAlignment} --maf2=${unmaskedAlignment} --out=${unmaskedQuality}

#Sampled alignments

${seedMultiplicityProbTable}: ${seedCounts} ${deBruijnGraph} ${globalSeq} ${seedProbabilitiesXX} ${seedProbabilitiesDir}
	${seedProbabilitiesXX} --seedMultiplicity --graphFile ${deBruijnGraph} --lastzFile ${seedCounts} --sequenceFile ${globalSeq} --statsFile ${seedMultiplicitySamplingStats} > ${seedMultiplicityProbTable}

${seedMultiplicityAlignment} ${seedMultiplicityAlignmentStats}: ${seedMultiplicityProbTable} ${lastzXX} ${trimmedSeq} ${alignmentsDir}
	${lastzXX} ${sampledLastzArgs} --seedSamplingProbabilities=${seedMultiplicityProbTable} --stats=${seedMultiplicityAlignmentStats} > ${seedMultiplicityAlignment}

${seedMultiplicityQuality} : ${seedMultiplicityAlignment} ${maskedAlignment}
	${mafComparatorXX} --maf1=${maskedAlignment} --maf2=${seedMultiplicityAlignment} --out=${seedMultiplicityQuality}


${pathMultiplicityProbTable}: ${seedCounts} ${deBruijnGraph} ${globalSeq} ${seedProbabilitiesXX} ${seedProbabilitiesDir}
	${seedProbabilitiesXX} --pathMultiplicity --graphFile ${deBruijnGraph} --lastzFile ${seedCounts} --sequenceFile ${globalSeq} --statsFile ${pathMultiplicitySamplingStats} > ${pathMultiplicityProbTable}

${pathMultiplicityAlignment} ${pathMultiplicityAlignmentStats}: ${pathMultiplicityProbTable} ${lastzXX} ${trimmedSeq} ${alignmentsDir}
	${lastzXX} ${sampledLastzArgs} --seedSamplingProbabilities=${pathMultiplicityProbTable} --stats=${pathMultiplicityAlignmentStats} > ${pathMultiplicityAlignment}

${pathMultiplicityQuality} : ${pathMultiplicityAlignment} ${maskedAlignment}
	${mafComparatorXX} --maf1=${maskedAlignment} --maf2=${pathMultiplicityAlignment} --out=${pathMultiplicityQuality}

seedProbabilities: ${pathMultiplicityProbTable} ${seedMultiplicityProbTable}

seedProbabilitiesClean:
	rm ${pathMultiplicityProbTable} ${seedMultiplicityProbTable}

alignments: ${maskedAlignment} ${unmaskedAlignment} ${seedMultiplicityAlignment} ${pathMultiplicityAlignment}

sampledAlignments: ${seedMultiplicityAlignment} ${pathMultiplicityAlignment}

quality: ${unmaskedQuality} ${seedMultiplicityQuality} ${pathMultiplicityQuality}


plots: ${pathMultiplicitySamplingStats}
	Rscript debruijnpaths.R

benchmarks: alignments quality plots



clean:
	rm -rf ${benchmarksDir}


