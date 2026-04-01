#!/bin/bash
# 用法：
#   bash merge_sum_by_KOID.sh file1.txt file2.txt ...
# 输出：
#   merged.txt

awk '
BEGIN {
    OFS = "\t"
}

# 读取每个文件
FNR == 1 {
    # 文件头第二列作为样本名
    sample[++h] = $2
    current_file = FILENAME
    next
}

{
    koid = $1
    value = $2 + 0
    sum[current_file, koid] += value   # 对同一文件中重复 KOID 求和
    all_koid[koid] = 1                 # 记录所有出现过的 KOID
}

END {
    # 打印表头
    printf "KOID"
    for (i = 1; i <= h; i++) {
        printf "\t%s", sample[i]
    }
    printf "\n"

    # 按 KOID 打印
    for (k in all_koid) {
        printf "%s", k
        for (i = 1; i <= h; i++) {
            file = ARGV[i]
            key = file SUBSEP k
            if (key in sum) {
                printf "\t%s", sum[key]
            } else {
                printf "\tNA"
            }
        }
        printf "\n"
    }
}' "$@" > merged.txt
