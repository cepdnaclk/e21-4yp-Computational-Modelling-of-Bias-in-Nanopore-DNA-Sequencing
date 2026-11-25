/*
 * Phase 1 pipeline: composition bias across basecalling, trimming, QC,
 * host mapping, classification and BLAST of failed reads.
 *
 * Run with ./run_pipeline.sh (it builds design.csv first). Each stage runs once
 * per unique upstream combination and its outputs are hard-linked into
 *   <outdir>/<key=value>/<key=value>/.../{metadata.json, results.tsv, data files}
 */

params.design = null
params.outdir = null
params.fast5 = []
params.threads = 16
params.gpu_device = 'cuda:all'
params.guppy_bin = 'guppy_basecaller'
params.dorado_bin = 'dorado'
params.dorado_models = ''
params.qc_constants = [maxlen: 50000]
params.host_map = [min_mapq: 20]
params.blast = [subsample: 2000, seed: 11, evalue: '1e-5', max_target_seqs: 5]


process INVENTORY {
    label 'basecall'
    publishDir "${params.outdir}/_shared/inventory", mode: 'link'

    input:
    path fast5, stageAs: 'fast5_in/*'

    output:
    path 'inventory.csv'
    path 'raw_reads.tsv.gz'

    script:
    """
    make_inventory.py --fast5 fast5_in/*.fast5 --outdir .
    """
}


process BASECALL_GUPPY {
    tag "${p}"
    label 'basecall'
    label 'gpu'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(p), val(spec)
    path fast5, stageAs: 'fast5_in/*'

    output:
    tuple val(p), path('reads.fastq.gz'), emit: reads
    path 'sequencing_summary.tsv.gz'
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${params.guppy_bin} -i fast5_in -s guppy_out -c ${spec} -x ${params.gpu_device} --disable_qscore_filtering --compress_fastq
    find guppy_out -name '*.fastq.gz' | sort | xargs cat > reads.fastq.gz
    gzip -c guppy_out/sequencing_summary.txt > sequencing_summary.tsv.gz
    EOF
    bash cmd.sh
    read_stats.py --reads reads.fastq.gz --outdir .
    write_metadata.py --stage basecall --path '${p}' --command-file cmd.sh \\
        --tool "guppy_basecaller=\$(${params.guppy_bin} --version | grep -o 'Version [^ ,]*' | head -1)" \\
        --param model_config=${spec} --param min_qscore=none --param device=${params.gpu_device} \\
        --param adapter_trimming=guppy_default --inputs ${fast5} \\
        --outputs reads.fastq.gz sequencing_summary.tsv.gz reads_table.tsv.gz results.tsv
    """
}


process BASECALL_DORADO {
    tag "${p}"
    label 'basecall'
    label 'gpu'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(p), val(spec)
    path fast5, stageAs: 'fast5_in/*'

    output:
    tuple val(p), path('reads.fastq.gz'), emit: reads
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    // Dorado 0.9.x reads FAST5 directly. (pod5 convert fails on these files: the
    // pip vbz plugin cannot decompress them, so POD5 is not used.)
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${params.dorado_bin} basecaller ${params.dorado_models}/${spec} fast5_in --device ${params.gpu_device} --min-qscore 0 --no-trim --emit-fastq | gzip -1 > reads.fastq.gz
    EOF
    bash cmd.sh
    read_stats.py --reads reads.fastq.gz --outdir .
    write_metadata.py --stage basecall --path '${p}' --command-file cmd.sh \\
        --tool "dorado=\$(${params.dorado_bin} --version 2>&1 | tail -1)" \\
        --param model=${spec} --param min_qscore=0 --param device=${params.gpu_device} \\
        --param adapter_trimming=none --param input_format=fast5 --inputs ${fast5} \\
        --outputs reads.fastq.gz reads_table.tsv.gz results.tsv
    """
}


process IMPORT_BASECALLS {
    tag "${p}"
    label 'basecall'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(p), val(dirs)

    output:
    tuple val(p), path('reads.fastq.gz'), emit: reads
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    def dir_list = dirs.tokenize(';').collect { d -> "'${d}'" }.join(' ')
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    find ${dir_list} -maxdepth 1 -name '*.fastq*' | sort | xargs zcat -f | gzip -1 > reads.fastq.gz
    EOF
    bash cmd.sh
    read_stats.py --reads reads.fastq.gz --outdir .
    write_metadata.py --stage basecall --path '${p}' --command-file cmd.sh \\
        --param source=imported_existing_basecalls --inputs ${dir_list} \\
        --outputs reads.fastq.gz reads_table.tsv.gz results.tsv
    """
}


process TRIM {
    tag "${p}"
    label 'qc'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(trim), path(reads, stageAs: 'input.fastq.gz')

    output:
    tuple val(p), path('reads.fastq.gz'), emit: reads
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    def cmd = trim == 'porechop_abi' ? "porechop_abi -abi -i input.fastq.gz -o reads.fastq.gz -t ${task.cpus}"
            : trim == 'porechop' ? "porechop -i input.fastq.gz -o reads.fastq.gz -t ${task.cpus}"
            : trim == 'none' ? "ln \"\$(readlink -f input.fastq.gz)\" reads.fastq.gz"
            : error("Unknown trim level: ${trim}")
    def ver = trim == 'none' ? '' : "--tool \"${trim}=\$(${trim} --version 2>&1 | tail -1)\""
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${cmd}
    EOF
    bash cmd.sh
    read_stats.py --reads reads.fastq.gz --outdir .
    write_metadata.py --stage trim --path '${p}' --parent '${parent}' --command-file cmd.sh ${ver} \\
        --inputs input.fastq.gz --outputs reads.fastq.gz reads_table.tsv.gz results.tsv
    """
}


process QC_SPLIT {
    tag "${p}"
    label 'qc'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(minq), val(minlen), path(reads, stageAs: 'input.fastq.gz')

    output:
    tuple val(p), path('passed.fastq.gz'), path('failed.fastq.gz'), emit: reads
    path 'artefact.fastq.gz'
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    qc_split.py --reads input.fastq.gz --minq ${minq} --minlen ${minlen} --maxlen ${params.qc_constants.maxlen} --outdir .
    EOF
    bash cmd.sh
    write_metadata.py --stage qc --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --param maxlen=${params.qc_constants.maxlen} \\
        --inputs input.fastq.gz \\
        --outputs passed.fastq.gz failed.fastq.gz artefact.fastq.gz reads_table.tsv.gz results.tsv
    """
}


process HOST_MAP {
    tag "${p}"
    label 'align'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(aligner), val(ref_name), path(ref), path(passed), path(failed)

    output:
    path 'passed.bam'
    path 'failed.bam'
    path 'host_reads.tsv.gz'
    path 'results.tsv', emit: results
    path 'metadata.json'

    script:
    if (aligner != 'minimap2') {
        error("Unknown aligner: ${aligner}")
    }
    def t = task.cpus
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    minimap2 -ax map-ont -t ${t} ${ref} ${passed} | samtools sort -@ 4 -o passed.bam
    minimap2 -ax map-ont -t ${t} ${ref} ${failed} | samtools sort -@ 4 -o failed.bam
    samtools index passed.bam
    samtools index failed.bam
    EOF
    bash cmd.sh
    host_map_summary.py --bam passed=passed.bam --bam failed=failed.bam --min-mapq ${params.host_map.min_mapq}
    write_metadata.py --stage host_map --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --tool "minimap2=\$(minimap2 --version)" --tool "samtools=\$(samtools --version | head -1)" \\
        --param host_ref=${ref_name} --param min_mapq=${params.host_map.min_mapq} \\
        --inputs ${ref} ${passed} ${failed} --outputs passed.bam failed.bam host_reads.tsv.gz results.tsv
    """
}


process CLASSIFY {
    tag "${p}"
    label 'classify'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(classifier), val(db), val(db_path), path(passed), path(failed)

    output:
    path '*.kreport'
    path '*.classified.gz'
    path 'top_taxa.tsv'
    path 'results.tsv', emit: results
    path 'metadata.json'

    script:
    def run_passed = classifyCmd(classifier, db_path, task.cpus, 'passed')
    def run_failed = classifyCmd(classifier, db_path, task.cpus, 'failed')
    def ver = classifier == 'kraken2' ? "kraken2 --version | head -1" : "centrifuge --version | head -1"
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${run_passed}
    ${run_failed}
    EOF
    bash cmd.sh
    kraken_summary.py --report passed=passed.kreport --report failed=failed.kreport
    write_metadata.py --stage classify --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --tool "${classifier}=\$(${ver})" --param db=${db} --param db_path=${db_path} \\
        --inputs ${passed} ${failed} --outputs passed.kreport failed.kreport top_taxa.tsv results.tsv
    """
}


process BLAST {
    tag "${p}"
    label 'blast'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(tool), val(db), val(db_prefix), val(db_fasta), path(passed), path(failed)

    output:
    path '*.sub.fa'
    path '*.hits'
    path 'best_hits.tsv'
    path 'results.tsv', emit: results
    path 'metadata.json'

    script:
    def b = params.blast
    def search_passed = searchCmd(tool, db_prefix, db_fasta, b, task.cpus, 'passed')
    def search_failed = searchCmd(tool, db_prefix, db_fasta, b, task.cpus, 'failed')
    def fmt = tool == 'minimap2' ? 'paf' : 'blast6'
    def ver = tool == 'minimap2' ? 'minimap2 --version' : 'blastn -version | head -1'
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    seqkit shuffle -s ${b.seed} passed.fastq.gz -o passed.shuf.fq.gz && seqkit head -n ${b.subsample} passed.shuf.fq.gz | seqkit fq2fa > passed.sub.fa && rm passed.shuf.fq.gz
    seqkit shuffle -s ${b.seed} failed.fastq.gz -o failed.shuf.fq.gz && seqkit head -n ${b.subsample} failed.shuf.fq.gz | seqkit fq2fa > failed.sub.fa && rm failed.shuf.fq.gz
    ${search_passed}
    ${search_failed}
    EOF
    bash cmd.sh
    blast_summary.py --format ${fmt} --hits passed=passed.hits --hits failed=failed.hits \\
        --queries passed=passed.sub.fa --queries failed=failed.sub.fa
    write_metadata.py --stage blast --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --tool "${tool}=\$(${ver})" --tool "seqkit=\$(seqkit version)" \\
        --param blast_db=${db} --param subsample=${b.subsample} --param seed=${b.seed} --param evalue=${b.evalue} \\
        --inputs ${passed} ${failed} --outputs passed.sub.fa failed.sub.fa passed.hits failed.hits best_hits.tsv results.tsv
    """
}


workflow {
    if (!params.design || !params.outdir) {
        error("Run via run_pipeline.sh (needs --design and --outdir)")
    }
    def fast5 = channel.fromPath(params.fast5, checkIfExists: true).collect()
    def rows = channel.fromPath(params.design).splitCsv(header: true)

    INVENTORY(fast5)

    // ---- basecall: one task per basecaller:model -----------------------------
    def bc = rows.map { r -> [r.bc_path, r.bc_type, r.model_spec, r.import_fastq] }.unique()
    BASECALL_GUPPY(bc.filter { r -> r[1] == 'guppy' }.map { r -> [r[0], r[2]] }, fast5)
    BASECALL_DORADO(bc.filter { r -> r[1] == 'dorado' }.map { r -> [r[0], r[2]] }, fast5)
    IMPORT_BASECALLS(bc.filter { r -> r[1] == 'import' }.map { r -> [r[0], r[3]] })
    def basecalled = BASECALL_GUPPY.out.reads.mix(BASECALL_DORADO.out.reads, IMPORT_BASECALLS.out.reads)

    // ---- trim: one task per basecall x trim ----------------------------------
    def trim_in = rows.map { r -> [r.bc_path, r.trim_path, r.trim] }.unique().combine(basecalled, by: 0)
    TRIM(trim_in)

    // ---- qc: one task per trim x minq x minlen -------------------------------
    def qc_in = rows.map { r -> [r.trim_path, r.qc_path, r.minq, r.minlen] }.unique().combine(TRIM.out.reads, by: 0)
    QC_SPLIT(qc_in)
    def qc_reads = QC_SPLIT.out.reads

    // ---- analyses: one task per leaf combination -----------------------------
    def host_in = rows.filter { r -> r.analysis == 'host_map' }
        .map { r -> [r.qc_path, r.path, r.aligner, r.host_ref, file(r.host_ref_fasta)] }
        .combine(qc_reads, by: 0)
    HOST_MAP(host_in)

    def cls_in = rows.filter { r -> r.analysis == 'classify' }
        .map { r -> [r.qc_path, r.path, r.classifier, r.db, r.db_path] }
        .combine(qc_reads, by: 0)
    CLASSIFY(cls_in)

    def blast_in = rows.filter { r -> r.analysis == 'blast' }
        .map { r -> [r.qc_path, r.path, r.blast_tool, r.blast_db, r.blast_db_prefix, r.blast_db_fasta] }
        .combine(qc_reads, by: 0)
    BLAST(blast_in)
}


// ---- command builders (one line per read set, so metadata shows the literal command)
def classifyCmd(classifier, db_path, t, s) {
    if (classifier == 'kraken2') {
        return "kraken2 --db ${db_path} --threads ${t} --gzip-compressed --report ${s}.kreport --output - ${s}.fastq.gz | gzip -1 > ${s}.classified.gz"
    }
    if (classifier == 'centrifuge') {
        return "centrifuge -x ${db_path} -U ${s}.fastq.gz -p ${t} -S ${s}.classified --report-file ${s}.centrifuge_report.tsv && centrifuge-kreport -x ${db_path} ${s}.classified > ${s}.kreport && gzip -1 ${s}.classified"
    }
    error("Unknown classifier: ${classifier}")
}

def searchCmd(tool, db_prefix, db_fasta, b, t, s) {
    if (tool == 'minimap2') {
        return "minimap2 -c -x map-ont -t ${t} ${db_fasta} ${s}.sub.fa > ${s}.hits"
    }
    if (tool in ['megablast', 'dc-megablast', 'blastn']) {
        return "blastn -task ${tool} -query ${s}.sub.fa -db ${db_prefix} -outfmt '6 qseqid sseqid pident length evalue bitscore' -evalue ${b.evalue} -max_target_seqs ${b.max_target_seqs} -num_threads ${t} -out ${s}.hits"
    }
    error("Unknown blast_tool: ${tool}")
}
