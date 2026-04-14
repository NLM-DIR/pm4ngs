import os

from pm4ngs.jupyterngsplugin.files.fastq.fastqc import parse_fastqc_zip
from pm4ngs.jupyterngsplugin.markdown.utils import find_file_print_link_size

from typing import List, Dict, Tuple


def fastqc_table(
        sample_list: List[str],
        samples_path: str,
        fastqc_path: str
) -> Tuple[Dict[str, dict], str]:
    """
    Generate a markdown table summarizing FastQC report results for a list of samples.

    Args:
        sample_list: List of sample prefix names.
        samples_path: Path to original fastq files.
        fastqc_path: Path to FastQC output folder.

    Returns:
        Tuple of (samples_data, table_string):
            samples_data: Dict with parsed FastQC data for each sample.
            table_string: Markdown-formatted table as a string.
    """
    samples_data: Dict[str, dict] = {}
    table_lines = [
        '| Sample | Fastq | FastQC<br>Report | No of Reads<br>in fastq | Seq<br> Len | %GC | Poor<br>Quality | Fail<br>Tests | Adapters |',
        '| --- | --- | --- |--- | --- | --- | --- | --- | --- |'
    ]
    for sample in sample_list:
        row = [sample]
        row.append(find_file_print_link_size(samples_path, sample, '.fastq.gz', 'MB', ' --- '))
        # Skip Jupyter checkpoint files (e.g., *-checkpoint.html)
        if not sample.endswith('-checkpoint'):
            row.append(find_file_print_link_size(fastqc_path, sample, '.html', 'MB', ' --- '))
        else:
            row.append(' --- ')
        fastqc_zip_path = os.path.relpath(os.path.join(fastqc_path, sample + '_fastqc.zip'))
        if os.path.exists(fastqc_zip_path) and os.path.getsize(fastqc_zip_path) != 0:
            tests, tot_seq, poor_quality, seq_len, gc_content, adapter = parse_fastqc_zip(fastqc_zip_path)
            samples_data[sample] = {'fastqc': {
                'tests': tests,
                'tot_seq': tot_seq,
                'poor_quality': poor_quality,
                'seq_len': seq_len,
                'gc_content': gc_content,
                'adapter': adapter
            }}
            row.append("{:,}".format(tot_seq))
            row.append(seq_len)
            row.append(gc_content)
            row.append(str(poor_quality))
            fail_tests = '<br>'.join([t for t in tests if tests[t] == 'FAIL'])
            row.append(fail_tests)
            row.append(adapter)
        else:
            row.extend([' --- '] * 6)
        table_lines.append('| ' + ' | '.join(row) + ' |')
    return samples_data, '\n'.join(table_lines)


def fastqc_trimmomatic_table(
        samples_data: Dict[str, dict],
        sample_list: List[str],
        trimmomatic_path: str
) -> Tuple[Dict[str, dict], str]:
    """
    Generate a markdown table summarizing FastQC report results for trimmed samples.

    Args:
        samples_data: Dict returned from fastqc_table with extracted data.
        sample_list: List of sample prefix names.
        trimmomatic_path: Folder with trimmed samples and FastQC reports.

    Returns:
        Tuple of (samples_data, table_string):
            samples_data: Dict updated with trimmed FastQC data for each sample.
            table_string: Markdown-formatted table as a string.
    """
    table_lines = [
        '| Sample | Trimmed<br>Fastq | FastQC<br>Report | No of Reads<br>in fastq | Removed<br>Reads | Seq<br> Len | %GC | Poor<br>Quality | Fail<br>Tests | Adapters |',
        '| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |'
    ]
    for sample in sample_list:
        row = [sample]
        row.append(find_file_print_link_size(trimmomatic_path, sample, '.fastq.gz', 'MB', ' --- '))
        row.append(find_file_print_link_size(trimmomatic_path, sample, '.html', 'KB', ' --- '))
        # Find trimmed fastqc zip
        files = [
            f for ds, dr, files in os.walk(trimmomatic_path) for f in files
            if f.startswith(sample) and f.endswith('_fastqc.zip')
               and os.path.getsize(os.path.join(trimmomatic_path, f)) != 0
        ]
        if len(files) == 1:
            fastqc_zip_path = os.path.relpath(os.path.join(trimmomatic_path, files[0]))
            tests, tot_seq, poor_quality, seq_len, gc_content, adapter = parse_fastqc_zip(fastqc_zip_path)
            row.append("{:,}".format(tot_seq))
            removed_reads = samples_data[sample]['fastqc']['tot_seq'] - tot_seq \
                if sample in samples_data and 'fastqc' in samples_data[sample] else 0
            row.append("{:,}".format(removed_reads))
            row.append(seq_len)
            row.append(gc_content)
            row.append(str(poor_quality))
            fail_tests = '<br>'.join([t for t in tests if tests[t] == 'FAIL'])
            row.append(fail_tests)
            row.append(adapter)
            samples_data[sample]['trimmed'] = {
                'tests': tests,
                'tot_seq': tot_seq,
                'poor_quality': poor_quality,
                'seq_len': seq_len,
                'gc_content': gc_content,
                'adapter': adapter
            }
        else:
            row.extend([' --- '] * 7)
        table_lines.append('| ' + ' | '.join(row) + ' |')
    return samples_data, '\n'.join(table_lines)
