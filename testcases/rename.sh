#!/usr/bin/env bash

# 找出所有 in/ans 檔案、取出前兩位數，排序，大排到小
for f in $(ls | grep -E '^[0-9]+\.((in)|(ans))$' | sort -r); do

    n=$(echo "$f" | sed -E 's/^([0-9]+)\..*/\1/')
    ext=${f##*.}
    new=$((10#$n + 2))

    printf -v new_name "%02d.%s" "$new" "$ext"

    echo "mv \"$f\" \"$new_name\""
    mv "$f" "$new_name"
done
