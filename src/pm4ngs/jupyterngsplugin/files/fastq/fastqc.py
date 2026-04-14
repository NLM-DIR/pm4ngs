import os
import pandas
import zipfile
from io import StringIO



def parse_fastqc_zip(zip_file):
    """
    This function extract from a fastQC zip file the failing tests and data
    :param zip_file: FastQC zip file
    :return:
    """
    tests = {}
    tot_seq = 0
    poor_quality = 0
    seq_len = ''
    gc_content = ''
    adapter_text = []
    if os.path.exists(zip_file) and os.path.getsize(zip_file) != 0:
        with zipfile.ZipFile(zip_file) as myzip:
            prefix = next(
                (name.split('/')[0] for name in myzip.namelist() if name.endswith('summary.txt')),
                None
            )
            with myzip.open(prefix + '/summary.txt') as myfile:
                for line in myfile:
                    fields = line.strip().decode('utf8').split('\t')
                    tests[fields[1]] = fields[0]
            with myzip.open(prefix + '/fastqc_data.txt') as myfile:
                adapter_flag = False
                for line in myfile:
                    line = line.strip().decode('utf8')
                    if adapter_flag and not line.startswith('>>'):
                        adapter_text.append(line.strip())
                    elif adapter_flag and line.startswith('>>'):
                        adapter_flag = False
                    fields = line.split('\t')
                    if fields[0] == 'Total Sequences':
                        tot_seq = int(fields[1])
                    if fields[0] == 'Sequences flagged as poor quality':
                        poor_quality = int(fields[1])
                    if fields[0] == 'Sequence length':
                        seq_len = fields[1]
                    if fields[0] == '%GC':
                        gc_content = fields[1]
                    if fields[0] == '>>Adapter Content':
                        adapter_flag = True
    if len(adapter_text) > 0:
        adapter_text = '\n'.join(adapter_text)
        df = pandas.read_csv(StringIO(adapter_text), sep="\t")
        df['#Position'] = df['#Position'].astype(str)
        averages = df.max(numeric_only=True)
        filtered_sorted = averages[averages > 1].sort_values(ascending=False)
        result_df = filtered_sorted.to_frame(name="average").reset_index(names="column")
        if not result_df.empty:
            adapter_text = f"{result_df.iloc[0]['column']} ({result_df.iloc[0]['average']:.2f})"
        else:
            adapter_text = ''
    else:
        adapter_text = ''
    return tests, tot_seq, poor_quality, seq_len, gc_content, adapter_text
