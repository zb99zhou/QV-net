#!/bin/bash

# Configuration
export NV=42
export NC=3
export TOKEN_NUM=13600000
TEST_CMD="cargo test test_voting_with_verify --release -- --nocapture"
RUN_TIMES=1

# Initialize data structures
declare -A OPERATIONS_COUNTS
declare -A OPERATIONS_SUMS
declare -A OPERATIONS_UNITS

reset_stats() {
    OPERATIONS_SUMS=()
    OPERATIONS_COUNTS=()
    OPERATIONS_UNITS=()
}

# Function: Convert all times to microseconds for consistent calculation
convert_to_us() {
    local value=$1
    local unit=$2
    
    case $unit in
        "ms") echo "$value * 1000" | bc -l ;;
        "µs") echo "$value" | bc -l ;;
        "ns") echo "$value / 1000" | bc -l ;;
        *) echo "$value" | bc -l ;;  # default to µs
    esac
}

# Function: Execute test and record timings
run_test() {
    local output=$($TEST_CMD 2>&1)
    # echo $output
    
    # Process each timing line
    while IFS= read -r line; do
        if [[ $line =~ "Time elapsed in "(.*)": "([0-9.]+)([^[:space:]]*) ]]; then
            local operation="${BASH_REMATCH[1]}"
            local value="${BASH_REMATCH[2]}"
            local unit="${BASH_REMATCH[3]}"
            
            # Convert to microseconds
            local value_us=$(convert_to_us "$value" "$unit")
            
            # Store original unit (use the first encountered)
            if [[ -z "${OPERATIONS_UNITS[$operation]}" ]]; then
                OPERATIONS_UNITS["$operation"]="$unit"
            fi
            
            # Accumulate values
            OPERATIONS_SUMS["$operation"]=$(echo "${OPERATIONS_SUMS[$operation]:-0} + $value_us" | bc -l)
            OPERATIONS_COUNTS["$operation"]=$(( ${OPERATIONS_COUNTS[$operation]:-0} + 1 ))
        fi
    done <<< "$output"
}

# Function: Calculate and display averages
display_results() {
    # echo "──────────────────────────────────────────────────────"
    # echo " TEST RESULTS SUMMARY (${RUN_TIMES} runs)"
    # echo "──────────────────────────────────────────────────────"
    # printf "%-25s %-12s %-12s\n" "OPERATION" "AVERAGE" "UNIT"
    echo "──────────────────────────────────────────────────────"
    
    for operation in "${!OPERATIONS_SUMS[@]}"; do
        local count=${OPERATIONS_COUNTS[$operation]}
        local sum_us=${OPERATIONS_SUMS[$operation]}
        local avg_us=$(echo "scale=6; $sum_us / $count" | bc -l)
        local unit=${OPERATIONS_UNITS[$operation]}
        # echo "$operation, $count, $sum_us, $avg_us, $unit"
        
        # Convert back to original unit for display
        case $unit in
            "ms") avg_display=$(echo "scale=6; $avg_us / 1000" | bc -l) ;;
            "ns") avg_display=$(echo "scale=6; $avg_us * 1000" | bc -l) ;;
            *) avg_display=$avg_us ;;
        esac
        
        printf "%-25s %-12.6f %-12s\n" "$operation" "$avg_display" "$unit"
    done | sort
    echo "──────────────────────────────────────────────────────"
}

# Main execution
cd ..
reset_stats
export NV=2
export NC=4
export TOKEN_NUM=100
echo "[Figure 3a - 3d]"
for ((i=1; i<=20; i++)); do
    reset_stats
    export NC=$i
    echo "Running with NC=$NC: Starting benchmark (${RUN_TIMES} iterations)..."
    for ((j=1; j<=$RUN_TIMES; j++)); do
        run_test
    done
    display_results
done



