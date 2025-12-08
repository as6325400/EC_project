#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $0 <path-to-solution-binary> [case-id]"
    exit 1
}

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TARGET_INPUT="$1"
CASE_FILTER="${2:-}"

if [[ "${TARGET_INPUT}" = /* ]]; then
    TARGET="${TARGET_INPUT}"
else
    TARGET="${ROOT_DIR}/${TARGET_INPUT}"
fi

if [ ! -x "${TARGET}" ]; then
    echo "Error: solution binary '${TARGET}' not found or not executable." >&2
    exit 1
fi

CHECKER="${ROOT_DIR}/bin/checker"
if [ ! -x "${CHECKER}" ]; then
    echo "Building checker..."
    make -C "${ROOT_DIR}" "${CHECKER}"
fi

TESTCASE_DIR="${ROOT_DIR}/testcases"
if [ ! -d "${TESTCASE_DIR}" ]; then
    echo "Error: testcases directory '${TESTCASE_DIR}' not found." >&2
    exit 1
fi

if [ -n "${CASE_FILTER}" ]; then
    candidate="${TESTCASE_DIR}/${CASE_FILTER}.in"
    if [ ! -f "${candidate}" ]; then
        echo "Error: testcase '${CASE_FILTER}' not found under '${TESTCASE_DIR}'." >&2
        exit 1
    fi
    test_inputs=("${candidate}")
else
    shopt -s nullglob
    test_inputs=("${TESTCASE_DIR}"/*.in)
    shopt -u nullglob
    if [ ${#test_inputs[@]} -eq 0 ]; then
        echo "Error: no .in files found under '${TESTCASE_DIR}'." >&2
        exit 1
    fi
fi

tmp_dir="$(mktemp -d)"
trap 'rm -rf "${tmp_dir}"' EXIT

passed=0
failed=0
total_score=0

printf "%-8s | %-12s | %-12s | %-10s\n" "Case" "Result" "Score" "Time(s)"
printf "%-8s-+-%-12s-+-%-12s-+-%-10s\n" "--------" "------------" "------------" "----------"

for input in "${test_inputs[@]}"; do
    case_name="$(basename "${input}")"
    case_name="${case_name%.in}"
    answer="${TESTCASE_DIR}/${case_name}.ans"

    result="AC"
    score_display="-"
    time_display="-"
    log_to_print=""

    if [ ! -f "${answer}" ]; then
        result="MissingAns"
        failed=$((failed + 1))
        printf "%-8s | %-12s | %-12s | %-10s\n" "${case_name}" "${result}" "${score_display}" "${time_display}"
        continue
    fi

    output_file="${tmp_dir}/${case_name}.out"
    log_file="${tmp_dir}/${case_name}.log"
    time_file="${tmp_dir}/${case_name}.time"

    if /usr/bin/env time -p -o "${time_file}" "${TARGET}" < "${input}" > "${output_file}"; then
        exec_status=0
    else
        exec_status=$?
    fi

    if [ -f "${time_file}" ]; then
        time_display="$(awk '/^real /{print $2}' "${time_file}" 2>/dev/null | head -n1)"
        if [ -z "${time_display}" ]; then
            time_display="-"
        fi
    fi

    if [ "${exec_status}" -ne 0 ]; then
        result="RE(${exec_status})"
        failed=$((failed + 1))
        printf "%-8s | %-12s | %-12s | %-10s\n" "${case_name}" "${result}" "${score_display}" "${time_display}"
        continue
    fi

    set +e
    "${CHECKER}" "${input}" "${answer}" < "${output_file}" > "${log_file}" 2>&1
    checker_status=$?
    set -e

    case "${checker_status}" in
        42)
            score_line="$(grep -m1 'OPT_SCORE=' "${log_file}" || true)"
            if [ -n "${score_line}" ]; then
                score_value="${score_line#*=}"
                if [[ "${score_value}" =~ ^-?[0-9]+$ ]]; then
                    score_display="${score_value}"
                    total_score=$((total_score + score_value))
                else
                    score_display="${score_line}"
                fi
            fi
            passed=$((passed + 1))
            ;;
        43)
            result="WA"
            log_to_print="${log_file}"
            failed=$((failed + 1))
            ;;
        *)
            result="CheckerErr(${checker_status})"
            log_to_print="${log_file}"
            failed=$((failed + 1))
            ;;
    esac

    printf "%-8s | %-12s | %-12s | %-10s\n" "${case_name}" "${result}" "${score_display}" "${time_display}"
    if [ -n "${log_to_print}" ]; then
        sed 's/^/    /' "${log_to_print}"
    fi
done

total=$((passed + failed))
printf "%-8s-+-%-12s-+-%-12s-+-%-10s\n" "--------" "------------" "------------" "----------"
printf "%-8s | %-12s | %-12s | %-10s\n" "Total" "${passed}/${total} passed" "${total_score}" "-"

if [ "${failed}" -ne 0 ]; then
    exit 1
fi
