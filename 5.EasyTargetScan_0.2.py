#!/usr/bin/python

# Easy TargetScan.py

from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import os
import re
import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
import argparse

# ========== 新增错误处理模块 ==========
import sys
import traceback

def log_error(error_msg):
    """记录错误信息到日志文件"""
    with open("EasyTargetScan_error.log", "a") as f:
        f.write(f"[ERROR] {error_msg}\n")
        traceback.print_exc(file=f)

# ========== 参数解析部分保持不变 ==========
parser = argparse.ArgumentParser(description='Allows to use TargetScan from FASTA files')
parser.add_argument("-utr", type=str, help="a FASTA file containing sequences to scan", required=True)
parser.add_argument("-mir", type=str, help="a FASTA file containing mature miRNA sequences", required=True)
parser.add_argument("-output", type=str, help="prefix of the output file", default="test_EasyTargetScan")
parser.add_argument("-Taxon", type=int, help="mouse by default:10090", default=10090)
args = parser.parse_args()

# ========== 核心改进部分 ==========
def main():
    try:
        # ========== 参数初始化 ==========
        utr_Database = args.utr
        mirna = args.mir
        seed_size = 7
        taxon = args.Taxon
        output_file = args.output
        path_to_targetscan = 'targetscan_70.pl'

        # ========== miRNA种子处理（使用安全写入方式）==========
        with open("miR_seeds_temp.txt", 'w') as f1:
            slen = seed_size + 1
            for miR in SeqIO.parse(mirna, "fasta"):
                try:
                    mir_seed = Seq(str(miR.seq)[1:slen])
                    f1.write(f"{miR.id}\t{str(mir_seed)}\t{taxon}\n")
                except IndexError:
                    print(f"警告: miRNA序列过短 ({miR.id})，已跳过")
                    continue

        # ========== UTR处理循环（关键改进部分）==========
        for utr in SeqIO.parse(utr_Database, "fasta"):
            try:
                # 改进点1: 过滤短序列
                if len(utr.seq) < 20:
                    print(f"跳过短UTR序列: {utr.id} (长度={len(utr.seq)})")
                    continue

                # 改进点2: 清理非法字符
                clean_utr_id = re.sub(r'[^a-zA-Z0-9_-]', '_', utr.id)
                output_file_utr = f"{clean_utr_id}_{output_file}"

                # 改进点3: 安全写入UTR临时文件
                with open("UTRs_temp.txt", 'w') as f2:
                    f2.write(f"{utr.id}\t{taxon}\t{utr.seq}\n")

                # 改进点4: 增强Perl执行错误处理
                cmd = f"perl {path_to_targetscan} miR_seeds_temp.txt UTRs_temp.txt {output_file_utr}"
                print(f"\n正在处理UTR: {utr.id}")
                print(f"执行命令: {cmd}")
                
                p = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True  # 确保文本模式输出
                )
                stdout, _ = p.communicate()

                # 检查执行状态
                if p.returncode != 0:
                    print(f"错误: Perl执行失败 (返回码={p.returncode})")
                    print("错误日志:")
                    print(stdout)
                    continue

                # 改进点5: 验证输出文件
                if not os.path.exists(output_file_utr):
                    print(f"严重错误: 未生成输出文件 {output_file_utr}")
                    continue

                # ========== 结果可视化处理 ==========
                features = []
                with open(output_file_utr) as fp:
                    next(fp)  # 跳过标题行
                    for line in fp:
                        try:
                            content = line.strip().split("\t")
                            features.append(
                                GraphicFeature(
                                    start=int(content[3]),
                                    end=int(content[4]),
                                    strand=+1,
                                    color=get_color(content[8]),
                                    label=re.sub(r'mmu-', '', content[1])
                                )
                            )
                        except Exception as e:
                            print(f"解析行错误: {line.strip()} | 错误: {str(e)}")
                            continue

                # ========== 生成可视化图表 ==========
                record = GraphicRecord(
                    sequence_length=len(str(utr.seq)),
                    features=features
                )
                record.plot(figure_width=12)
                plt.title(f'{utr.id} sequence')
                
                # 创建图例
                legend_handles = [
                    mpatches.Patch(color="#00ff99", label='6mer'),
                    mpatches.Patch(color="#9999ff", label='7mer-1a'),
                    mpatches.Patch(color="#ff66cc", label='7mer-m8'),
                    mpatches.Patch(color="#ff0000", label='8mer-1a')
                ]
                plt.legend(handles=legend_handles)
                
                # 自动保存图片
                plt.savefig(f"{output_file_utr}.png", dpi=300)
                plt.close()  # 防止内存泄漏

            except Exception as e:
                log_error(f"处理UTR {utr.id} 时发生错误: {str(e)}")
                continue

    finally:
        # ========== 清理临时文件 ==========
        temp_files = ["miR_seeds_temp.txt", "UTRs_temp.txt"]
        for f in temp_files:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except Exception as e:
                    print(f"清理临时文件失败: {f} | 错误: {str(e)}")

def get_color(seed_match_type):
    return {
        '6mer': "#00ff99",
        '7mer-1a': "#9999ff",
        '7mer-m8': "#ff66cc",
        '8mer-1a': "#ff0000"
    }.get(seed_match_type, "#ccccff")

if __name__ == "__main__":
    main()
