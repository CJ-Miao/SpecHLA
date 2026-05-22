"""
HLA Typing from diploid assemblies using BLAST
基于 annoHLA.pl 的 BLAST 比对和评分逻辑

输入：两个单倍型组装序列（fasta 格式）
输出：HLA 判型结果（使用 BLAST 评分，更接近 SpecHLA 主流程）

用法:
python3 typing_from_assembly_blast.py -n sample -1 hap1.fa -2 hap2.fa -o output --gene A
"""

import os
import sys
import re
import argparse
import subprocess
import tempfile
from collections import defaultdict
from spechla_paths import get_db_dir
import pysam


def extract_regions(input_fasta, gene, hap_index, output_fasta):
    """
    根据 annoHLA.pl 的逻辑提取特定区域
    不同基因有不同的多样性区域
    """
    # 定义各基因的区域 (start-end, 基于 1 的坐标)
    regions = {
        'A': [(100, 3300)],
        'B': [(150, 4000)],
        'C': [(400, 3500)],
        'DPA1': [(400, 700), (4200, 5500)],
        'DPB1': [(300, 600), (5000, 5300), (9200, 10600)],
        'DQA1': [(500, 900), (4600, 6100)],
        'DQB1': [(500, 2400), (5200, 7300)],
        'DRB1': [(100, 11000)],
    }

    if gene not in regions:
        # 如果没有定义区域，使用完整序列
        with open(input_fasta, 'r') as fin:
            with open(output_fasta, 'w') as fout:
                for line in fin:
                    fout.write(line)
        return

    try:
        # 使用 pysam 读取 fasta 文件
        with pysam.FastaFile(input_fasta) as fa:
            contigs = fa.references
            if not contigs:
                # 如果没有 contig 名，使用默认名
                contig = 'HLA_{}_{}'.format(gene, hap_index - 1)
            else:
                contig = contigs[0]

            with open(output_fasta, 'w') as fout:
                for start, end in regions[gene]:
                    try:
                        seq = fa.fetch(contig, start - 1, end)  # pysam 使用 0-based 坐标
                        fout.write(f">{contig}:{start}-{end}\n{seq}\n")
                    except ValueError:
                        # 如果区域超出范围，使用可用区域
                        seq_len = fa.get_reference_length(contig)
                        actual_start = max(0, start - 1)
                        actual_end = min(seq_len, end)
                        if actual_start < actual_end:
                            seq = fa.fetch(contig, actual_start, actual_end)
                            fout.write(f">{contig}:{start}-{end}\n{seq}\n")
    except Exception as e:
        print(f"    Warning: Could not extract regions: {e}")
        # 如果提取失败，使用完整序列
        with open(input_fasta, 'r') as fin:
            with open(output_fasta, 'w') as fout:
                for line in fin:
                    fout.write(line)


def get_IMGT_version():
    """获取 IMGT/HLA 数据库版本"""
    g_file = f"{get_db_dir()}/HLA/hla_nom_g.txt"
    version_info = "N/A"
    try:
        for line in open(g_file):
            if "# version:" in line:
                version_info = line.strip()
                break
    except:
        pass
    return version_info


def run_blast(query_fasta, db_path, out_file, threads=4):
    """运行 BLAST 比对"""
    # 检查数据库是否存在
    if not os.path.exists(db_path + ".nsq"):
        print(f"Error: BLAST database not found: {db_path}")
        return False

    cmd = f"blastn -query {query_fasta} -db {db_path} -out {out_file} -outfmt 7 -num_threads {threads} -max_target_seqs 6000"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"BLAST error: {result.stderr}")
        return False
    return True


def write_blast_summary(blast_file, hap_index, summary_file, gene):
    """
    写入 BLAST 汇总文件，格式与 SpecHLA 的 hla.blast.summary.txt 一致
    格式：HLA_{gene}_{hap}  allele_name  mismatch+gap  total_len  score
    """
    allele_scores = defaultdict(lambda: {'target_len': 0, 'mismatch_gap': 0})

    try:
        with open(blast_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                hla_allele = fields[1]
                if not hla_allele.startswith(f"{gene}*"):
                    continue

                try:
                    target_len = int(fields[3])
                    mismatch = int(fields[4])
                    gap_open = int(fields[5])

                    allele_scores[hla_allele]['target_len'] += target_len
                    allele_scores[hla_allele]['mismatch_gap'] += mismatch + gap_open
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        return

    # 写入汇总文件
    with open(summary_file, 'a') as f:
        for allele, data in sorted(allele_scores.items(), key=lambda x: x[0]):
            if data['target_len'] > 0:
                score = 100 * (1 - data['mismatch_gap'] / data['target_len'])
                # 格式：HLA_A_1  A*01:01:01:01  mismatch+gap  total_len  score
                f.write(f"HLA_{gene}_{hap_index}\t{allele}\t{data['mismatch_gap']}\t{data['target_len']}\t{score:.2f}\n")


def parse_blast_result(blast_file, gene):
    """
    解析 BLAST 结果，计算每个等位基因的得分
    逻辑来自 annoHLA.pl whole_blast:
    - hash11{hla} += target_length (比对总长度)
    - hash12{hla} += mismatch + gap (错配 + gap)
    - score = 100 * (1 - hash12/hash11)
    实际上就是比对一致性 (identity)
    """
    allele_scores = defaultdict(lambda: {'target_len': 0, 'mismatch_gap': 0, 'hits': 0})

    try:
        with open(blast_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                # outfmt 7 字段: query, subject, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
                # 我们使用 outfmt 6 的字段: 1=subject, 3=pident, 4=length, 5=mismatch, 6=gapopen
                hla_allele = fields[1]  # subject (等位基因名)

                # 只处理目标基因
                if not hla_allele.startswith(f"{gene}*"):
                    continue

                try:
                    target_len = int(fields[3])    # 比对长度
                    mismatch = int(fields[4])      # 错配数
                    gap_open = int(fields[5])      # gap 数

                    # annoHLA.pl: hash12 = mismatch + gap
                    allele_scores[hla_allele]['target_len'] += target_len
                    allele_scores[hla_allele]['mismatch_gap'] += mismatch + gap_open
                    allele_scores[hla_allele]['hits'] += 1
                except (ValueError, IndexError) as e:
                    continue
    except FileNotFoundError:
        print(f"Warning: BLAST result file not found: {blast_file}")
        return []

    # 计算得分并排序
    results = []
    for allele, data in allele_scores.items():
        if data['hits'] == 0 or data['target_len'] == 0:
            continue

        # score = 100 * (1 - mismatch_gap/target_len)
        # 这就是比对一致性 (identity)
        score = 100 * (1 - data['mismatch_gap'] / data['target_len'])

        results.append({
            'allele': allele,
            'score': score,
            'target_len': data['target_len'],
            'mismatch_gap': data['mismatch_gap'],
            'hits': data['hits']
        })

    # 按得分排序
    results.sort(key=lambda x: x['score'], reverse=True)
    return results


def load_population_freq(db_dir):
    """
    加载人群频率数据库
    格式：Gene*C*08:01  Caucasian_freq  Black_freq  Asian_freq
    """
    freq_file = f"{db_dir}/HLA/HLA_FREQ_HLA_I_II.txt"
    freq_dict = {}
    try:
        with open(freq_file, 'r') as f:
            next(f)  # 跳过表头
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    allele = parts[0]  # 如 C*08:01
                    freq_dict[allele] = {
                        'caucasian': float(parts[1]),
                        'black': float(parts[2]),
                        'asian': float(parts[3])
                    }
    except FileNotFoundError:
        print(f"Warning: Population frequency file not found: {freq_file}")
    return freq_dict


def select_best_allele(results, min_score=50, pop_freq_dict=None, pop='Unknown'):
    """
    选择最佳等位基因
    逻辑来自 annoHLA.pl whole_blast:
    - 得分 >= min_score 的候选
    - 过滤掉没有人群频率信息的等位基因（除非 pop='nonuse'）
    - 选择得分最高的（如果有多个得分相同，选择字母顺序最靠后的）
    """
    if not results:
        return None

    # 过滤掉没有频率信息的等位基因（除非 pop='nonuse'）
    filtered_results = []
    for r in results:
        allele = r['allele']
        # 提取前两位（如 A*01:01:01:01 -> A*01:01）
        parts = allele.split(':')
        two_field = f"{parts[0]}:{parts[1]}" if len(parts) >= 2 else allele

        if pop != 'nonuse':
            # 检查是否有频率信息
            if pop_freq_dict and two_field in pop_freq_dict:
                freq_info = pop_freq_dict[two_field]
                # 检查是否有任何人群频率 > 0
                if freq_info['caucasian'] > 0 or freq_info['black'] > 0 or freq_info['asian'] > 0:
                    filtered_results.append(r)
            elif pop_freq_dict is None:
                # 如果没有加载频率数据库，不过滤
                filtered_results.append(r)
        else:
            filtered_results.append(r)

    if not filtered_results:
        # 如果所有等位基因都被过滤了，使用原始结果
        filtered_results = results

    # 按得分排序，然后按字母顺序排序（字母顺序靠后的优先）
    filtered_results.sort(key=lambda x: (x['score'], x['allele']), reverse=True)

    best = filtered_results[0]
    if best['score'] < min_score:
        return None

    # 返回所有得分相同的最佳等位基因
    best_score = best['score']
    best_alleles = [r for r in filtered_results if abs(r['score'] - best_score) < 0.01]

    return best_alleles


def main():
    parser = argparse.ArgumentParser(
        description="HLA Typing from diploid assemblies using BLAST",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    required.add_argument("-n", type=str, help="Sample ID", metavar="")
    required.add_argument("-1", type=str, help="Assembly file of haplotype 1", metavar="")
    required.add_argument("-2", type=str, help="Assembly file of haplotype 2", metavar="")
    required.add_argument("-o", type=str, help="Output directory", default="./output")

    optional.add_argument("--gene", type=str, help="Specify gene(s) to type, e.g., 'A' or 'A,B,C'", default=None)
    optional.add_argument("-j", type=int, help="Number of threads", default=4)
    optional.add_argument("--db", type=str, help="Database directory", default=get_db_dir())
    optional.add_argument("--min-score", type=float, help="Minimum score threshold", default=50.0)
    optional.add_argument("--pop", type=str, help="Population: Asian, Black, Caucasian, Unknown, nonuse", default='Unknown')

    args = vars(parser.parse_args())

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    # 确定要分析的基因
    all_genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    if args.get('gene'):
        gene_list = [g.strip() for g in args['gene'].split(',')]
        gene_list = [g for g in gene_list if g in all_genes]
        if not gene_list:
            print("Error: No valid genes specified")
            sys.exit(1)
    else:
        gene_list = all_genes

    print(f"Analyzing genes: {', '.join(gene_list)}")

    # 创建输出目录
    result_path = args['o']
    os.makedirs(result_path, exist_ok=True)

    sample = args['n']
    hap1_file = args['1']
    hap2_file = args['2']
    db_dir = args['db']
    threads = args['j']
    min_score = args['min_score']
    pop = args['pop']

    version_info = get_IMGT_version()
    print(f"Database version: {version_info}")
    print(f"Population: {pop}")

    # 加载人群频率数据库
    pop_freq_dict = load_population_freq(db_dir)

    # 创建 BLAST 汇总文件（与 SpecHLA 格式一致）
    blast_summary_file = f"{result_path}/hla.blast.summary.txt"
    # 如果文件已存在则先删除（避免追加重复内容）
    if os.path.exists(blast_summary_file):
        os.remove(blast_summary_file)

    # 临时目录
    tmpdir = tempfile.mkdtemp(dir=result_path)

    results = {1: {}, 2: {}}  # hap1_results, hap2_results

    for hap_index, hap_file in [(1, hap1_file), (2, hap2_file)]:
        print(f"\nProcessing haplotype {hap_index}...")

        for gene in gene_list:
            print(f"  Typing HLA-{gene}...")

            # 准备 query 文件（提取特定区域，与 annoHLA.pl 一致）
            query_file = f"{tmpdir}/hap{hap_index}_{gene}.fa"

            # 提取特定区域（与 annoHLA.pl 一致）
            extract_regions(hap_file, gene, hap_index, query_file)

            # BLAST 数据库路径
            # 对于 DRB1，当 pop != "nonuse" 时使用 exon 数据库（与 annoHLA.pl 一致）
            if gene == "DRB1" and pop != "nonuse":
                db_path = f"{db_dir}/HLA/whole/HLA_DRB1.exon"
                if not os.path.exists(db_path + ".nsq"):
                    db_path = f"{db_dir}/HLA/whole/HLA_DRB1.exon.fa"
                    if not os.path.exists(db_path + ".nsq"):
                        db_path = f"{db_dir}/HLA/whole/HLA_DRB1"
            else:
                db_path = f"{db_dir}/HLA/whole/HLA_{gene}"

            if not os.path.exists(db_path + ".nsq"):
                # 尝试备用路径
                db_path = f"{db_dir}/HLA/whole/HLA_{gene}.fa"
                if not os.path.exists(db_path + ".nsq"):
                    db_path = f"{db_dir}/HLA/whole/HLA_{gene}.fasta"
                    if not os.path.exists(db_path + ".nsq"):
                        # 尝试使用 hla_gen.fasta
                        db_path = f"{db_dir}/HLA/whole/hla_gen"
                        if not os.path.exists(db_path + ".nsq"):
                            print(f"    Error: No BLAST database found for HLA_{gene}")
                            continue

            # 运行 BLAST
            blast_out = f"{tmpdir}/hap{hap_index}_{gene}.blast"
            if run_blast(query_file, db_path.replace('.fa', '').replace('.fasta', ''), blast_out, threads):
                # 写入 BLAST 汇总文件（与 SpecHLA 格式一致）
                write_blast_summary(blast_out, hap_index, blast_summary_file, gene)

                # 解析结果
                blast_results = parse_blast_result(blast_out, gene)

                if blast_results:
                    best = select_best_allele(blast_results, min_score, pop_freq_dict, pop)
                    if best:
                        results[hap_index][gene] = best
                        print(f"    Best allele: {best[0]['allele']} (score: {best[0]['score']:.2f})")
                    else:
                        print(f"    No allele passed score threshold")
                        results[hap_index][gene] = []
                else:
                    print(f"    No BLAST hits found")
                    results[hap_index][gene] = []
            else:
                print(f"    BLAST failed")
                results[hap_index][gene] = []

    # 保留临时目录中的 BLAST 中间文件（移动到结果目录）
    blast_dir = f"{result_path}/blast_tmp"
    os.makedirs(blast_dir, exist_ok=True)
    for gene in gene_list:
        for hap in [1, 2]:
            src = f"{tmpdir}/hap{hap}_{gene}.blast"
            dst = f"{blast_dir}/hap{hap}_{gene}.blast"
            if os.path.exists(src):
                subprocess.run(f"mv {src} {dst}", shell=True)
    # 删除临时目录的其他文件
    subprocess.run(f"rm -rf {tmpdir}", shell=True)

    # 输出结果
    result_file = f"{result_path}/hla.typing.results.txt"
    with open(result_file, 'w') as f:
        f.write("Locus\tChromosome\tAllele\tIdentity\n")
        for gene in gene_list:
            for hap in [1, 2]:
                if gene in results[hap] and results[hap][gene]:
                    best = results[hap][gene][0]
                    # Convert score (0-100) to identity (0-1)
                    identity = best['score'] / 100.0
                    f.write(f"{gene}\t{hap}\t{best['allele']}\t{identity:.4f}\n")
                else:
                    f.write(f"{gene}\t{hap}\tN/A\tN/A\n")

    print(f"\nTyping results saved to {result_file}")

    # 打印摘要
    print("\n" + "="*50)
    print("TYPING SUMMARY")
    print("="*50)
    for gene in gene_list:
        h1_allele = results[1].get(gene, [{}])[0].get('allele', 'N/A') if results[1].get(gene) else 'N/A'
        h2_allele = results[2].get(gene, [{}])[0].get('allele', 'N/A') if results[2].get(gene) else 'N/A'
        h1_identity = results[1].get(gene, [{}])[0].get('score', 'N/A') if results[1].get(gene) else 'N/A'
        h2_identity = results[2].get(gene, [{}])[0].get('score', 'N/A') if results[2].get(gene) else 'N/A'
        if isinstance(h1_identity, float):
            h1_identity = f"{h1_identity/100:.4f}"
        if isinstance(h2_identity, float):
            h2_identity = f"{h2_identity/100:.4f}"
        print(f"HLA-{gene}: {h1_allele} (identity: {h1_identity}) / {h2_allele} (identity: {h2_identity})")


if __name__ == "__main__":
    main()
