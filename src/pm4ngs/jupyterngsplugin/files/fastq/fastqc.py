import os
import pandas
import zipfile
from io import StringIO


def _read_summary(myzip, prefix: str) -> dict:
    tests: dict[str, str] = {}
    with myzip.open(prefix + '/summary.txt') as myfile:
        for raw_line in myfile:
            fields = raw_line.strip().decode('utf8').split('\t')
            tests[fields[1]] = fields[0]
    return tests


def _parse_fastqc_data(myzip, prefix: str) -> tuple[int, int, str, str, list[str]]:
    total_sequences = 0
    poor_quality_sequences = 0
    sequence_length = ''
    gc_content = ''
    adapter_lines: list[str] = []
    adapter_section = False

    with myzip.open(prefix + '/fastqc_data.txt') as myfile:
        for raw_line in myfile:
            line = raw_line.strip().decode('utf8')

            if adapter_section:
                if line.startswith('>>'):
                    adapter_section = False
                else:
                    adapter_lines.append(line)
                    continue

            fields = line.split('\t')
            key = fields[0]

            if key == 'Total Sequences':
                total_sequences = int(fields[1])
            elif key == 'Sequences flagged as poor quality':
                poor_quality_sequences = int(fields[1])
            elif key == 'Sequence length':
                sequence_length = fields[1]
            elif key == '%GC':
                gc_content = fields[1]
            elif key == '>>Adapter Content':
                adapter_section = True

    return total_sequences, poor_quality_sequences, sequence_length, gc_content, adapter_lines


def _process_adapter_content(adapter_lines: list[str]) -> str:
    if not adapter_lines:
        return ''

    adapter_text = '\n'.join(adapter_lines)
    df = pandas.read_csv(StringIO(adapter_text), sep='\t')
    df['#Position'] = df['#Position'].astype(str)

    max_values = df.max(numeric_only=True)
    filtered = max_values[max_values > 1].sort_values(ascending=False)

    if filtered.empty:
        return ''

    top = filtered.iloc[0]
    column = filtered.index[0]
    return f"{column} ({top:.2f})"


def parse_fastqc_zip(zip_file: str) -> tuple[dict, int, int, str, str, str]:
    """
    Extract key metrics and failing tests from a FastQC ZIP file.

    Args:
        zip_file (str): Path to the FastQC ZIP file.

    Returns:
        tuple:
            - dict: Test status (test name -> PASS/WARN/FAIL)
            - int: Total sequences
            - int: Poor quality sequences
            - str: Sequence length
            - str: GC content
            - str: Dominant adapter contamination summary
    """
    if not (os.path.exists(zip_file) and os.path.getsize(zip_file) > 0):
        return {}, 0, 0, '', '', ''

    with zipfile.ZipFile(zip_file) as myzip:
        prefix = next(
            (name.split('/')[0] for name in myzip.namelist() if name.endswith('summary.txt')),
            None
        )

        if prefix is None:
            return {}, 0, 0, '', '', ''

        tests = _read_summary(myzip, prefix)

        (
            total_sequences,
            poor_quality_sequences,
            sequence_length,
            gc_content,
            adapter_lines,
        ) = _parse_fastqc_data(myzip, prefix)

    adapter_summary = _process_adapter_content(adapter_lines)

    return (
        tests,
        total_sequences,
        poor_quality_sequences,
        sequence_length,
        gc_content,
        adapter_summary,
    )
