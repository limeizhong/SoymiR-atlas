#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import gzip
import sys

def read_even_lines(filename):
    even_lines = []
    with gzip.open(filename, 'rt', encoding='utf-8') as file:
        lines = file.readlines()
        even_lines = [line.strip() for line in lines[1::2]]  # 从索引1开始，步长为2，读取偶数行
    return even_lines

def count_characters_and_length(char_list):
    char_info = {}  # 用字典来存储每个字符的频数和长度
    for char in char_list:
        # 保留 19-24 nt 的 sRNA 序列，并记录每条序列的出现次数。
        if 19 <= len(char) <= 24:
            if char not in char_info:
                char_info[char] = {"frequency": 1, "length": len(char)}
            else:
                # 如果字符已经在字典中，更新频数
                char_info[char]["frequency"] += 1

    return char_info

if __name__ == '__main__': 
    filename=sys.argv[1]
    output_filename = filename.replace('.fa.gz', '.filter.fa')
    even_lines=read_even_lines(filename)
    result = count_characters_and_length(even_lines)
    id=0
    with open(output_filename, 'wt', encoding='utf-8') as output_file:
        for char, info in result.items():
            output_file.write(f">read{id}_x{info['frequency']}\n")
            output_file.write(f"{char}\n")
            id=id+1
