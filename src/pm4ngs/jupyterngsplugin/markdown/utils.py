import os

from pm4ngs.jupyterngsplugin.image.pdftobase64 import pdftobase64
from pm4ngs.jupyterngsplugin.image.pngtobase64 import imagetobase64


def hide_code_str():
    """
    Print an HTML code to hide python code in Jupyter notebooks
    :return: String
    """
    return '''
    <script>
        code_show=true;
        function code_toggle() {
            if (code_show){
                $('div.input').hide();
            } else {
                $('div.input').show();
            }
            code_show = !code_show
        }
        $( document ).ready(code_toggle);
    </script>
    The raw code for this IPython notebook is by default hidden for easier reading.
    To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>.
    '''


def get_link_size(filename, size_unit="MB", ifempty=""):
    """
    Print an HTML link for the file using the size as text
    :param filename: file to be linked
    :param size_unit: Size units to display KB, MG or GB.
           It will display bytes for any other option
    :param ifempty: Use this if the output is empty
    :return: HTML string
    """
    str_msg = ifempty
    if os.path.exists(filename):
        size = os.path.getsize(filename)
        if size != 0:
            if size_unit == "KB":
                size = size / 1024.0
            elif size_unit == "MB":
                size = size / (1024.0 ** 2)
            elif size_unit == "GB":
                size = size / (1024.0 ** 3)
            str_msg = ' <a href="'
            str_msg += filename.replace(' ', '%20')
            str_msg += '" target="_blank">'
            str_msg += "{0:.2f}".format(size)
            str_msg += size_unit
            str_msg += '</a>'
    return str_msg


def get_link_image(filename, width, height, ifempty=''):
    """
    Print an HTML link for the file using ta reduce image
    :param filename: file to be linked
    :param width: Image width
    :param height: Image height
    :param ifempty: Use this if the output is empty
    :return: HTML string
    """
    str_msg = ifempty
    if os.path.exists(filename) and os.path.getsize(filename) != 0:
        str_msg = ' <a href="'
        str_msg += filename.replace(' ', '%20')
        str_msg += '" target="_blank">'
        if filename.endswith('.pdf'):
            str_msg += '<img src="data:image/png;base64,' + pdftobase64(filename, width, height)
        else:
            str_msg += '<img src="data:image/png;base64,' + imagetobase64(filename, width, height)
        str_msg += '">'
    return str_msg


def get_link_name(filename, ifempty=""):
    """
    Print an HTML link for the file using the basename as text
    :param filename: file to be linked
    :param ifempty: Use this if the output is empty
    :return: HTML string
    """
    str_msg = ifempty
    if os.path.exists(filename) and os.path.getsize(filename) != 0:
        str_msg = ' <a href="'
        str_msg += filename.replace(' ', '%20')
        str_msg += '" target="_blank">'
        str_msg += os.path.basename(filename)
        str_msg += '</a>'
    return str_msg


def get_link_text(filename, text, ifempty=""):
    """
    Print an HTML link for the file using the basename as text
    :param filename: file to be linked
    :param text: Text to use in the link
    :param ifempty: Use this if the output is empty
    :return: HTML string
    """
    str_msg = ifempty
    if os.path.exists(filename) and os.path.getsize(filename) != 0:
        str_msg = ' <a href="'
        str_msg += filename.replace(' ', '%20')
        str_msg += '" target="_blank">'
        str_msg += '{0}'.format(text)
        str_msg += '</a>'
    return str_msg


def find_file_print_link_size(folder_path, prefix, sufix, size_unit, ifempty):
    """
    Find file in a folder and print link with size
    :param folder_path: Folder containing the file
    :param prefix: File name prefix
    :param sufix: File name sufix
    :param size_unit: Size units to display KB, MG or GB.
           It will display bytes for any other option
    :param ifempty: Use this if the output is empty
    :return:
    """
    files = [
        f for ds, dr, files in os.walk(folder_path)
        for f in files
        if f.startswith(prefix) and f.endswith(sufix)
        and os.path.exists(os.path.join(folder_path, f))
        and os.path.getsize(os.path.join(folder_path, f)) != 0
    ]
    if len(files) == 1:
        f = os.path.relpath(os.path.join(folder_path, files[0]))
        return get_link_size(f, size_unit, ifempty)
    return ifempty


def find_file_print_link_name(folder_path, prefix, sufix, size_unit, ifempty):
    """
    Find file in a folder and print link with name
    :param folder_path: Folder containing the file
    :param prefix: File name prefix
    :param sufix: File name sufix
    :param size_unit: Size units to display KB, MG or GB.
           It will display bytes for any other option
    :param ifempty: Use this if the output is empty
    :return:
    """
    files = [f for ds, dr, files in os.walk(folder_path)
             for f in files if f.startswith(prefix) and f.endswith(sufix)
             and os.path.getsize(os.path.join(folder_path, f)) != 0]
    if len(files) == 1:
        f = os.path.relpath(os.path.join(folder_path, files[0]))
        return get_link_name(f, ifempty)
    return ifempty


def info_table(notebook, result_dir):
    str_msg = '#### Info\n\n'
    str_msg += '| Type | Link |\n'
    str_msg += '| --- | --- |\n'
    str_msg += '| Notebook | <a href="' \
               + notebook.replace(' ', '%20') \
               + '.ipynb" target="_blank">' \
               + notebook + '</a> | \n'
    str_msg += '| Results | <a href="' \
               + os.path.relpath(result_dir).replace(' ', '%20') \
               + '" target="_blank">' \
               + os.path.basename(result_dir) + '</a> |\n\n'

    return str_msg
