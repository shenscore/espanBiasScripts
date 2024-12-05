#!/bin/bash

# 检查参数数量
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <regions.bed> <replication_origins.bed> <positive_strand.bw> <negative_strand.bw> <outputPrefix>"
    exit 1
fi

# 输入参数
REGIONS=$1
ORIGINS=$2
POS_BW=$3
NEG_BW=$4
RESULTS=$5
echo "pos file is $POS_BW"
echo "neg file is $NEG_BW"

# 输出中间文件及最终结果文件名
CLOSEST="closest_results.bed"
FILTERED="filtered_regions.bed"
POS_SIGNAL="positive_signal.tab"
NEG_SIGNAL="negative_signal.tab"
POS_BED="positive_regions.bed"
NEG_BED="negative_regions.bed"
POS_SUM="positive_sum.txt"
NEG_SUM="negative_sum.txt"
COMBINED="combined_signals.bed"
STRAND_ASSIGNMENT="strand_assignment.bed"

# Step 1: 计算最近复制起始位点及其距离
echo "Step 1: Calculating closest replication origins..."
bedtools closest -a "$REGIONS" -b "$ORIGINS" -D ref > "$CLOSEST"

# Step 2: 过滤不符合条件的区间
echo "Step 2: Filtering regions based on distance..."
awk '$NF <= 100000 && $NF != 0' "$CLOSEST" > "$FILTERED"

# Step 3: 为 bigWigAverageOverBed 准备标准的 4 列 BED 文件
echo "Preparing 4-column BED files for signal extraction..."
awk '{print $1, $2, $3, "region_"NR}' OFS="\t" "$FILTERED" > "$POS_BED"
cp "$POS_BED" "$NEG_BED"

# Step 4: 提取正负链信号值，并取第 4 列 sum 值
echo "Step 4: Extracting signal values from bigwig files..."
bigWigAverageOverBed "$POS_BW" "$POS_BED" "$POS_SIGNAL"
bigWigAverageOverBed "$NEG_BW" "$NEG_BED" "$NEG_SIGNAL"

# 提取第 4 列 sum 值，合并到过滤的区间文件
awk '{print ($4 < 0 ? -$4 : $4)}' "$POS_SIGNAL" > "$POS_SUM" # 取绝对值
awk '{print ($4 < 0 ? -$4 : $4)}' "$NEG_SIGNAL" > "$NEG_SUM"
paste "$FILTERED" "$POS_SUM" "$NEG_SUM" > "$COMBINED"

# Step 5: 根据上下游关系标记 leading 和 lagging 链, 这里的distance位置是repli origin相对于region的，所以要反过来
echo "Step 5: Assigning leading and lagging strands..."
awk -vOFS='\t' '{
    if ($(NF-2) < 0) {
        print $0, $(NF-1), $(NF); # 正链为 leading
    } else {
        print $0, $(NF), $(NF-1); # 负链为 leading
    }
}' "$COMBINED" > "$STRAND_ASSIGNMENT"

# Step 6: 计算 espan_bias
echo "Step 6: Calculating espan_bias..."

awk -vOFS='\t' '{
    leading = $(NF-1);
    lagging = $NF;
    if (leading + lagging == 0) {
        espan_bias = "NA"; # 避免除零错误，标记为 NA
    } else {
        espan_bias = (leading - lagging) / (leading + lagging);
    }
    print $0, espan_bias;
}' "$STRAND_ASSIGNMENT" > "$RESULTS"

# 清理中间文件
echo "Cleaning up intermediate files..."
rm -f "$CLOSEST" "$FILTERED" "$POS_SIGNAL" "$NEG_SIGNAL" "$POS_BED" "$NEG_BED" "$POS_SUM" "$NEG_SUM" "$COMBINED" "$STRAND_ASSIGNMENT"

# 输出结果
echo "Analysis complete. Results saved to $RESULTS."

