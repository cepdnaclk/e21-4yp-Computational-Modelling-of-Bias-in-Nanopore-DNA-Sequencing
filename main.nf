/*
 * Phase 1 pipeline: composition bias across basecalling, trimming, QC,
 * host mapping, classification and read search (BLAST and alternatives).
 *
 * Run with ./run_pipeline.sh (it builds design.csv first). Each stage runs once
 * per unique upstream combination and its outputs are hard-linked into
 *   <outdir>/<key=value>/<key=value>/.../{metadata.json, results.tsv, data files}
 * Tools that live outside a stage's own conda environment are added to PATH
 * from params.conda_base (set by run_pipeline.sh).
 */

params.design = null
params.outdir = null
params.fast5 = []
params.threads = 16
params.gpu_device = 'cuda:all'
params.bonito_batchsize = 64   // Bonito's default batch runs out of memory on a 15 GB A16
params.conda_base = "${System.getenv('HOME')}/anaconda3"
params.taxonomy = ''
params.qc_constants = [maxlen: 50000]
params.host_map = [mapq: '0,20,60']
params.trim = [lsk110_top: 'CCTGTACTTCGTTCAGTTACGTATTGCT', lsk110_bottom: 'GCAATACGTAACTGAACGAAGTACAGG',
               cutadapt_error: 0.2, cutadapt_overlap: 12, headcrop: 30, tailcrop: 10, scan_window: 150]
params.blast = [subsample: 2000, seed: 11, evalue: '1e-5', max_target_seqs: 5]
params.host_taxid = '3118'   // genus Ulva: the host clade when grouping taxonomy-based hits


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


process BASECALL {
    tag "${p}"
    label 'basecall'
    label 'gpu'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(p), val(bc_type), val(bc_bin), val(spec)
    path fast5, stageAs: 'fast5_in/*'

    output:
    tuple val(p), path('reads.fastq.gz'), emit: reads
    path 'sequencing_summary.tsv.gz', optional: true
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    def cmd = basecallCmd(bc_type, bc_bin, spec, params.gpu_device, params.conda_base, params.bonito_batchsize)
    def ver = bc_type == 'guppy' ? "${bc_bin} --version | grep -o 'Version [^ ,]*' | head -1"
            : bc_type == 'dorado' ? "${bc_bin} --version 2>&1 | tail -1"
            : "${bc_bin} --version 2>&1 | grep -v -i warn | tail -1"
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${ind(cmd)}
    EOF
    bash cmd.sh
    read_stats.py --reads reads.fastq.gz --outdir .
    write_metadata.py --stage basecall --path '${p}' --command-file cmd.sh \\
        --tool "${bc_type}=\$(${ver})" --param model=${spec} --param min_qscore=none \\
        --param device=${params.gpu_device} --param adapter_trimming=none_or_basecaller_default --inputs ${fast5} \\
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


process FASTQC {
    // The lab protocol quality-checked raw basecalls with FastQC; same report here.
    tag "${p}"
    label 'qc'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(p), path(reads, stageAs: 'reads.fastq.gz')

    output:
    path 'fastqc/*'

    script:
    """
    mkdir -p fastqc
    fastqc --threads ${task.cpus} --outdir fastqc --quiet reads.fastq.gz
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
    def cmd = trimCmd(trim, params.trim, task.cpus)
    def tool = trim == 'cutadapt_lsk110' ? 'cutadapt' : trim == 'headcrop' ? 'chopper' : trim == 'porechop_lsk110' ? 'porechop_abi' : trim
    def ver = trim == 'none' ? '' : "--tool \"${tool}=\$(${tool} --version 2>&1 | tail -1)\""
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${ind(cmd)}
    EOF
    bash cmd.sh
    read_stats.py --reads reads.fastq.gz --outdir .
    adapter_scan.py --reads reads.fastq.gz --top ${params.trim.lsk110_top} --bottom ${params.trim.lsk110_bottom} \\
        --window ${params.trim.scan_window} --max-error ${params.trim.cutadapt_error}
    write_metadata.py --stage trim --path '${p}' --parent '${parent}' --command-file cmd.sh ${ver} \\
        --param lsk110_top=${params.trim.lsk110_top} --param lsk110_bottom=${params.trim.lsk110_bottom} \\
        --param adapter_scan="edlib, window ${params.trim.scan_window} nt, max error ${params.trim.cutadapt_error}" \\
        --inputs input.fastq.gz --outputs reads.fastq.gz reads_table.tsv.gz results.tsv
    """
}


process QC_SPLIT {
    tag "${p}"
    label 'qc'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(qfilter), val(minq), val(keep), val(minlen), path(reads, stageAs: 'input.fastq.gz')

    output:
    tuple val(p), path('passed.fastq.gz'), path('failed.fastq.gz'), emit: reads
    path 'artefact.fastq.gz'
    path 'reads_table.tsv.gz'
    path 'results.tsv'
    path 'metadata.json'

    script:
    def maxlen = params.qc_constants.maxlen
    def cmd = qfilter == 'filtlong'
        ? "seqkit seq -M ${maxlen} input.fastq.gz -o under_max.fastq.gz\n" +
          "filtlong --keep_percent ${keep} --min_length ${Math.max(1, minlen as int)} under_max.fastq.gz | awk 'NR % 4 == 1 {print substr(\$1, 2)}' > pass_ids.txt\n" +
          "qc_split.py --reads input.fastq.gz --pass-ids pass_ids.txt --minlen ${minlen} --maxlen ${maxlen} --outdir .\n" +
          "rm under_max.fastq.gz"
        : "qc_split.py --reads input.fastq.gz --qmetric ${qfilter} --minq ${minq} --minlen ${minlen} --maxlen ${maxlen} --outdir ."
    def ver = qfilter == 'filtlong' ? "--tool \"filtlong=\$(filtlong --version 2>&1 | tail -1)\"" : ''
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    ${ind(cmd)}
    EOF
    bash cmd.sh
    write_metadata.py --stage qc --path '${p}' --parent '${parent}' --command-file cmd.sh ${ver} \\
        --param maxlen=${maxlen} --inputs input.fastq.gz \\
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
    path '*.bam', optional: true
    path 'host_reads.tsv.gz'
    path 'results.tsv', emit: results
    path 'metadata.json'

    script:
    def t = task.cpus
    def setup = aligner == 'winnowmap'
        ? "meryl count k=15 output merylDB ${ref}\nmeryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt"
        : aligner == 'lra' ? "lra index -ONT ${ref}" : ''
    def summary = aligner == 'hostile'
        ? '--clean passed=passed.fastq.gz,passed.clean.fastq.gz --clean failed=failed.fastq.gz,failed.clean.fastq.gz'
        : '--bam passed=passed.bam --bam failed=failed.bam'
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    export PATH="\$PATH:${params.conda_base}/envs/env-hostile/bin"
    ${ind(setup)}
    ${alignCmd(aligner, ref, t, 'passed')}
    ${alignCmd(aligner, ref, t, 'failed')}
    EOF
    bash cmd.sh
    host_map_summary.py ${summary} --mapq ${params.host_map.mapq}
    write_metadata.py --stage host_map --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --tool "${aligner}=\$(bash -c 'export PATH=\$PATH:${params.conda_base}/envs/env-hostile/bin; ${alignerVersion(aligner)}' 2>&1 | head -1)" \\
        --param host_ref=${ref_name} --param mapq_cutoffs=${params.host_map.mapq} \\
        --inputs ${ref} ${passed} ${failed} --outputs host_reads.tsv.gz results.tsv
    """
}


process CLASSIFY {
    tag "${p}"
    label 'classify'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(classifier), val(refdb), val(db_path), val(conf), path(passed), path(failed)

    output:
    path '*.kreport'
    path '*.out.gz', optional: true
    path 'top_taxa.tsv'
    path 'results.tsv', emit: results
    path 'metadata.json'

    script:
    def tax = params.taxonomy
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    export PATH="\$PATH:${params.conda_base}/envs/env-taxa/bin:${params.conda_base}/envs/env-protein/bin"
    ${ind(classifyCmd(classifier, db_path, conf, task.cpus, 'passed', tax))}
    ${ind(classifyCmd(classifier, db_path, conf, task.cpus, 'failed', tax))}
    for f in *.out; do [ -e "\$f" ] && gzip -f "\$f"; done
    EOF
    bash cmd.sh
    kraken_summary.py --report passed=passed.kreport --report failed=failed.kreport
    write_metadata.py --stage classify --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --tool "${classifier}=\$(bash -c 'export PATH=\$PATH:${params.conda_base}/envs/env-taxa/bin:${params.conda_base}/envs/env-protein/bin; ${classifierVersion(classifier)}' 2>&1 | grep -v -i warn | head -1)" \\
        --param refdb=${refdb} --param db_path=${db_path} --param confidence=${conf ?: 'na'} --param taxonomy=${tax} \\
        --inputs ${passed} ${failed} --outputs passed.kreport failed.kreport top_taxa.tsv results.tsv
    """
}


process BLAST {
    tag "${p}"
    label 'blast'
    publishDir path: { "${params.outdir}/${p}" }, mode: 'link'

    input:
    tuple val(parent), val(p), val(tool), val(refdb), val(db_prefix), val(db_fasta), path(passed), path(failed)

    output:
    path '*.sub.fa'
    path '*.hits'
    path 'best_hits.tsv'
    path 'results.tsv', emit: results
    path 'metadata.json'

    script:
    def b = params.blast
    def fmt = tool == 'minimap2' ? 'paf' : 'blast6'
    def ver = tool == 'minimap2' ? 'minimap2 --version' : tool == 'diamond' ? 'diamond version'
            : tool == 'last' ? 'lastal --version' : tool == 'mmseqs2' ? 'mmseqs version' : 'blastn -version | head -1'
    """
    cat > cmd.sh <<'EOF'
    set -euo pipefail
    seqkit shuffle -s ${b.seed} passed.fastq.gz -o passed.shuf.fq.gz && seqkit head -n ${b.subsample} passed.shuf.fq.gz | seqkit fq2fa > passed.sub.fa && rm passed.shuf.fq.gz
    seqkit shuffle -s ${b.seed} failed.fastq.gz -o failed.shuf.fq.gz && seqkit head -n ${b.subsample} failed.shuf.fq.gz | seqkit fq2fa > failed.sub.fa && rm failed.shuf.fq.gz
    ${searchCmd(tool, db_prefix, db_fasta, b, task.cpus, 'passed')}
    ${searchCmd(tool, db_prefix, db_fasta, b, task.cpus, 'failed')}
    EOF
    bash cmd.sh
    blast_summary.py --format ${fmt} --hits passed=passed.hits --hits failed=failed.hits \\
        --queries passed=passed.sub.fa --queries failed=failed.sub.fa --taxonomy ${params.taxonomy} --host-taxid ${params.host_taxid}
    write_metadata.py --stage blast --path '${p}' --parent '${parent}' --command-file cmd.sh \\
        --tool "${tool}=\$(${ver} 2>&1 | head -1)" --tool "seqkit=\$(seqkit version)" \\
        --param refdb=${refdb} --param db_prefix=${db_prefix} --param subsample=${b.subsample} --param seed=${b.seed} --param evalue=${b.evalue} \\
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
    def bc = rows.map { r -> [r.bc_path, r.bc_type, r.bc_bin, r.model_spec, r.import_fastq] }.unique()
    BASECALL(bc.filter { r -> r[1] != 'import' }.map { r -> [r[0], r[1], r[2], r[3]] }, fast5)
    IMPORT_BASECALLS(bc.filter { r -> r[1] == 'import' }.map { r -> [r[0], r[4]] })
    def basecalled = BASECALL.out.reads.mix(IMPORT_BASECALLS.out.reads)
    FASTQC(basecalled)

    // ---- trim: one task per basecall x trim ----------------------------------
    def trim_in = rows.map { r -> [r.bc_path, r.trim_path, r.trim] }.unique().combine(basecalled, by: 0)
    TRIM(trim_in)

    // ---- qc: one task per trim x qfilter x (minq | keep) x minlen -------------------
    def qc_in = rows.map { r -> [r.trim_path, r.qc_path, r.qfilter, r.minq, r.keep, r.minlen] }.unique()
        .combine(TRIM.out.reads, by: 0)
    QC_SPLIT(qc_in)
    def qc_reads = QC_SPLIT.out.reads

    // ---- analyses: one task per leaf combination -----------------------------
    def host_in = rows.filter { r -> r.analysis == 'host_map' }
        .map { r -> [r.qc_path, r.path, r.aligner, r.host_ref, file(r.host_ref_fasta)] }
        .combine(qc_reads, by: 0)
    HOST_MAP(host_in)

    def cls_in = rows.filter { r -> r.analysis == 'classify' }
        .map { r -> [r.qc_path, r.path, r.classifier, r.refdb, r.db_path, r.confidence] }
        .combine(qc_reads, by: 0)
    CLASSIFY(cls_in)

    def blast_in = rows.filter { r -> r.analysis == 'blast' }
        .map { r -> [r.qc_path, r.path, r.blast_tool, r.refdb, r.blast_db_prefix, r.blast_db_fasta] }
        .combine(qc_reads, by: 0)
    BLAST(blast_in)
}


// Multi-line commands are indented to match the script block, so Nextflow's
// indentation stripping leaves the heredoc's EOF marker at column 0.
def ind(s) {
    return s.toString().replace('\n', '\n    ')
}

// ==== command builders (one literal command per read set, so metadata shows it) =====

def basecallCmd(type, bin, spec, device, condaBase, bonitoBatch) {
    if (type == 'guppy') {
        return "${bin} -i fast5_in -s guppy_out -c ${spec} -x ${device} --disable_qscore_filtering --compress_fastq\n" +
               "find guppy_out -name '*.fastq.gz' | sort | xargs cat > reads.fastq.gz\n" +
               "gzip -c guppy_out/sequencing_summary.txt > sequencing_summary.tsv.gz"
    }
    if (type == 'dorado') {
        return "${bin} basecaller ${spec} fast5_in --device ${device} --min-qscore 0 --no-trim --emit-fastq | gzip -1 > reads.fastq.gz"
    }
    if (type == 'bonito') {
        def dev = device == 'cuda:all' ? 'cuda:0' : device
        return "export PATH=${condaBase}/envs/env-bonito/bin:\$PATH\n" +
               "${bin} basecaller ${spec} fast5_in --device ${dev} --no-trim --batchsize ${bonitoBatch} | gzip -1 > reads.fastq.gz"
    }
    error("Unknown basecaller type: ${type}")
}

def trimCmd(trim, tp, t) {
    if (trim == 'none') {
        return "ln \"\$(readlink -f input.fastq.gz)\" reads.fastq.gz"
    }
    if (trim == 'porechop') {
        return "porechop -i input.fastq.gz -o reads.fastq.gz -t ${t}"
    }
    if (trim == 'porechop_abi') {
        return "porechop_abi -abi -i input.fastq.gz -o reads.fastq.gz -t ${t}"
    }
    if (trim == 'porechop_lsk110') {
        // Porechop's algorithm with ONLY the lab protocol's LSK110 ligation adapter (no built-in list, no inference)
        return "printf 'LSK110_LA\\n${tp.lsk110_top}\\n${tp.lsk110_bottom}\\n' > lsk110_adapters.txt\n" +
               "porechop_abi -cap lsk110_adapters.txt -ddb -i input.fastq.gz -o reads.fastq.gz -t ${t}"
    }
    if (trim == 'cutadapt_lsk110') {
        // LSK110 ligation adapter (lab protocol): top strand at the 5' end, its complement at the 3' end
        return "cutadapt -g ${tp.lsk110_top} -a ${tp.lsk110_bottom} -e ${tp.cutadapt_error} -O ${tp.cutadapt_overlap} " +
               "--times 2 -j ${t} -o reads.fastq.gz input.fastq.gz > cutadapt.log"
    }
    if (trim == 'headcrop') {
        return "zcat input.fastq.gz | chopper --headcrop ${tp.headcrop} --tailcrop ${tp.tailcrop} --threads ${t} | gzip -1 > reads.fastq.gz"
    }
    error("Unknown trim level: ${trim}")
}

def alignCmd(aligner, ref, t, s) {
    def empty = "if [ -z \"\$(zcat ${s}.fastq.gz | head -c 1)\" ]; then samtools dict ${ref} | samtools view -b -o ${s}.bam -; "
    def sort = "samtools sort -@ 4 -o ${s}.bam && samtools index ${s}.bam"
    if (aligner == 'minimap2') {
        return empty + "else minimap2 -ax map-ont -t ${t} ${ref} ${s}.fastq.gz | ${sort}; fi"
    }
    if (aligner == 'winnowmap') {
        return empty + "else winnowmap -W repetitive_k15.txt -ax map-ont -t ${t} ${ref} ${s}.fastq.gz | ${sort}; fi"
    }
    if (aligner == 'ngmlr') {
        return empty + "else zcat ${s}.fastq.gz > ${s}.fq && ngmlr -t ${t} -x ont -r ${ref} -q ${s}.fq | ${sort} && rm ${s}.fq; fi"
    }
    if (aligner == 'lra') {
        return empty + "else zcat ${s}.fastq.gz > ${s}.fq && lra align -ONT ${ref} ${s}.fq -t ${t} -p s | ${sort} && rm ${s}.fq; fi"
    }
    if (aligner == 'hostile') {
        return "if [ -z \"\$(zcat ${s}.fastq.gz | head -c 1)\" ]; then cp ${s}.fastq.gz ${s}.clean.fastq.gz; " +
               "else hostile clean --fastq1 ${s}.fastq.gz --aligner minimap2 --index ${ref} --airplane --threads ${t} " +
               "--output hostile_${s} > hostile_${s}.json && mv hostile_${s}/*.clean*.fastq.gz ${s}.clean.fastq.gz; fi"
    }
    error("Unknown aligner: ${aligner}")
}

def alignerVersion(aligner) {
    def v = [minimap2: 'minimap2 --version', winnowmap: 'winnowmap --version', ngmlr: 'ngmlr --version',
             lra: 'lra -v', hostile: 'hostile --version']
    return v[aligner] ?: 'echo unknown'
}

def classifyCmd(clf, db, conf, t, s, tax) {
    def n = "N=\$(zcat -f ${s}.fastq.gz | awk 'END {print int(NR / 4)}')"
    def kr = "reads2kreport.py --total \$N --taxonomy ${tax} --out ${s}.kreport"
    def body
    if (clf == 'kraken2') {
        body = "kraken2 --db ${db} --threads ${t} --confidence ${conf ?: 0} --gzip-compressed --report ${s}.kreport --output ${s}.out ${s}.fastq.gz"
    } else if (clf == 'krakenuniq') {
        body = "krakenuniq --db ${db} --threads ${t} --gzip-compressed --report-file ${s}.native_report --output ${s}.out ${s}.fastq.gz && " +
               "awk -F'\\t' '\$1 == \"C\" {print \$2 \"\\t\" \$3}' ${s}.out > ${s}.assign && ${kr} --assign ${s}.assign"
    } else if (clf == 'centrifuge') {
        body = "centrifuge -x ${db}/index -U ${s}.fastq.gz -p ${t} -S ${s}.out --report-file ${s}.native_report && " +
               "centrifuge-kreport -x ${db}/index ${s}.out > ${s}.kreport"
    } else if (clf == 'centrifuger') {
        body = "centrifuger -x ${db}/index -u ${s}.fastq.gz -t ${t} > ${s}.out && centrifuger-kreport -x ${db}/index ${s}.out > ${s}.kreport"
    } else if (clf == 'kaiju') {
        body = "kaiju -t ${tax}/nodes.dmp -f ${db}/kaiju_db.fmi -i ${s}.fastq.gz -z ${t} -o ${s}.out && " +
               "awk -F'\\t' '\$1 == \"C\" {print \$2 \"\\t\" \$3}' ${s}.out > ${s}.assign && ${kr} --assign ${s}.assign"
    } else if (clf == 'metabuli') {
        body = "metabuli classify --seq-mode 3 --threads ${t} --max-ram 32 ${s}.fastq.gz ${db} mb_${s} ${s} > mb_${s}.log && " +
               "cp mb_${s}/${s}_classifications.tsv ${s}.out && awk -F'\\t' '\$1 == 1 {print \$2 \"\\t\" \$3}' ${s}.out > ${s}.assign && " +
               "${kr} --assign ${s}.assign"
    } else if (clf == 'ganon') {
        // ganon exits 1 when no read classifies ("No valid report entries"); accept only that case.
        body = "{ ganon classify -d ${db}/db -s ${s}.fastq.gz -o ${s} -t ${t} --output-one > ${s}.ganon.log 2>&1 || grep -q '^0 sequences classified' ${s}.ganon.log; } && touch ${s}.one && cp ${s}.one ${s}.out && " +
               "awk -F'\\t' 'NR == FNR {m[\$1] = \$2; next} {print \$1 \"\\t\" ((\$2 in m) ? m[\$2] : \$2)}' ${db}/target2taxid.tsv ${s}.one > ${s}.assign && " +
               "${kr} --assign ${s}.assign"
    } else if (clf == 'taxor') {
        // taxor keeps the whole FASTQ header as the read name, so give it read IDs only.
        // 'taxor profile' (0.2.1) crashes in its EM step on this data, so each read gets the
        // taxid of its best reference by matching k-mers (QHASH_MATCH) from 'taxor search'.
        body = "seqkit seq -i ${s}.fastq.gz -o ${s}.ids.fastq.gz && " +
               "taxor search --index-file ${db}/taxor.hixf --query-file ${s}.ids.fastq.gz --output-file ${s}.out --threads ${t} > /dev/null && " +
               "awk -F'\t' 'NR > 1 && NF >= 10 && \$2 != \"-\" { if (!(\$1 in b) || \$8 > b[\$1]) { b[\$1] = \$8; tx[\$1] = \$4 } } " +
               "END { for (r in tx) print r \"\\t\" tx[r] }' ${s}.out > ${s}.assign && rm ${s}.ids.fastq.gz && ${kr} --assign ${s}.assign"
    } else if (clf == 'kmcp') {
        body = "kmcp search -d ${db}/kmcp.db ${s}.fastq.gz -o ${s}.search.tsv.gz -j ${t} 2> ${s}.kmcp.log && " +
               "kmcp profile -X ${tax} -T ${db}/taxid.map ${s}.search.tsv.gz -o ${s}.profile -B ${s} -j ${t} 2>> ${s}.kmcp.log && " +
               "{ zcat -f ${s}.binning* 2>/dev/null || true; } > ${s}.out && ${kr} --cami ${s}.out"
    } else if (clf == 'sylph') {
        body = "sylph profile ${db}/db.syldb ${s}.fastq.gz -t ${t} > ${s}.out && ${kr} --sylph ${s}.out --genome-map ${db}/genome2taxid.tsv"
    } else if (clf == 'diamond') {
        body = "diamond blastx -d ${db}/db.dmnd -q ${s}.fastq.gz --long-reads -p ${t} --outfmt 102 -o ${s}.out --quiet && " +
               "awk -F'\\t' '\$2 != 0 {print \$1 \"\\t\" \$2}' ${s}.out > ${s}.assign && ${kr} --assign ${s}.assign"
    } else if (clf == 'megan') {
        body = "diamond blastx -d ${db}/db.dmnd -q ${s}.fastq.gz --long-reads -p ${t} --outfmt 100 -o ${s}.daa --quiet && " +
               "daa-meganizer -i ${s}.daa --longReads -alg longReads -a2t ${db}/acc2taxid.tsv > ${s}.megan.log 2>&1 && " +
               "daa2info -i ${s}.daa -r2c Taxonomy > ${s}.out 2>> ${s}.megan.log && " +
               "awk -F'\\t' '\$2 != \"\" && \$2 != 0 {print \$1 \"\\t\" \$2}' ${s}.out > ${s}.assign && ${kr} --assign ${s}.assign"
    } else {
        error("Unknown classifier: ${clf}")
    }
    // '|| exit 1': inside if/else, a failed step of an && chain would otherwise not stop the script
    return "${n}\nif [ \"\$N\" -eq 0 ]; then printf '0.00\\t0\\t0\\tU\\t0\\tunclassified\\n' > ${s}.kreport; " +
           "else ${body} || { echo '${clf} failed on ${s}' >&2; exit 1; }; fi"
}

def classifierVersion(clf) {
    def v = [kraken2: 'kraken2 --version', krakenuniq: 'krakenuniq --version', centrifuge: 'centrifuge --version',
             centrifuger: 'centrifuger -v', kaiju: 'kaiju -h 2>&1 | grep -i kaiju', metabuli: 'metabuli version',
             ganon: 'ganon --version', taxor: 'taxor --version', kmcp: 'kmcp version', sylph: 'sylph --version',
             diamond: 'diamond version', megan: 'daa-meganizer -h 2>&1 | grep -i version']
    return v[clf] ?: 'echo unknown'
}

def searchCmd(tool, db_prefix, db_fasta, b, t, s) {
    def guard = "if [ ! -s ${s}.sub.fa ]; then : > ${s}.hits; else "
    if (tool == 'minimap2') {
        return guard + "minimap2 -c -x map-ont -t ${t} ${db_fasta} ${s}.sub.fa > ${s}.hits; fi"
    }
    if (tool in ['megablast', 'dc-megablast', 'blastn']) {
        return guard + "BLASTDB=\$(dirname ${db_prefix}) blastn -task ${tool} -query ${s}.sub.fa -db ${db_prefix} -outfmt '6 std staxids' -evalue ${b.evalue} " +
               "-max_target_seqs ${b.max_target_seqs} -num_threads ${t} -out ${s}.hits; fi"
    }
    if (tool == 'diamond') {
        return guard + "diamond blastx -d ${db_prefix}.dmnd -q ${s}.sub.fa --long-reads -p ${t} --outfmt 6 " +
               "-e ${b.evalue} -k ${b.max_target_seqs} -o ${s}.hits --quiet; fi"
    }
    if (tool == 'last') {
        return guard + "lastal -P ${t} -f BlastTab ${db_prefix}.lastdb ${s}.sub.fa > ${s}.hits; fi"
    }
    if (tool == 'mmseqs2') {
        return guard + "mmseqs easy-search ${s}.sub.fa ${db_prefix}.mmseqs ${s}.hits mmseqs_tmp_${s} --search-type 3 " +
               "--threads ${t} -e ${b.evalue} > mmseqs_${s}.log && rm -rf mmseqs_tmp_${s}; fi"
    }
    error("Unknown blast_tool: ${tool}")
}
