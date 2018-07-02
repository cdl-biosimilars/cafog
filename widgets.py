import os
from PyQt5.QtWidgets import QFileDialog


class FileTypes:
    """
    A class for managing file types in :class:`~PyQt5.QtWidgets.QFileDialog` .

    :cvar dict _default_file_types: default file types
    :ivar list types: actually used file types
    :ivar list filters: filter strings

    .. automethod:: __init__
    """

    _default_file_types = {
        "xls": ("xls", "Excel 97-2003"),
        "xlsx": ("xlsx", "Excel"),
        "csv": ("csv", "Comma-separated value"),
        "": ("", "all files")
    }

    def __init__(self, type_list=None):
        """
        Create a new :class:`FileTypes` object.

        :param type_list: list of (extension, description) tuples
                          or {ext: description} dict
                          or list of extensions, which are then filtered
                          from ``_default_file_types``
        :type type_list: list(tuple) or dict or list(str)
        :return: nothing
        :raise TypeError: if an invalid type was specified for ``ext_list``
        """

        if type_list is None:
            self.types = self._default_file_types
            return

        self.types = []
        try:  # dict
            for ext, description in sorted(type_list.items()):
                self.types.append((ext, description))
        except AttributeError:
            try:  # list of tuples
                for ext, description in type_list:
                    self.types.append((ext, description))
            except (TypeError, ValueError):
                try:  # list of strings
                    for ext in type_list:
                        try:
                            self.types.append(self._default_file_types[ext])
                        except KeyError:
                            pass
                except (TypeError, ValueError):
                    raise TypeError("Argument for 'ext_list' has wrong type")
        self.filters = ["{1} [{0}] (*.{0})".format(k, v)
                        for k, v in self.types]

    def get_filter(self):
        """
        Returns  a file filter suitable for the ``filter`` parameter of
        :class:`~PyQt5.QtWidgets.QFileDialog`.

        :return: a file filter
        :rtype: str
        """

        return ";;".join(self.filters)

    def get_type(self, filefilter):
        """
        Returns the extension and description associated with a file filter.

        :param str filefilter: filefilter as returned by the dialog
        :return: the extension and description of the file type
        :rtype: tuple(str, str)
        """

        return self.types[self.filters.index(filefilter)]


def get_filename(parent, kind="save", caption="", directory="",
                 file_types=None):
    """
    Get a filename by a :class:`QFileDialog`
    and automatically add extensions.

    :param QWidget parent: parent of the dialog
    :param str kind: ``"save"`` or ``"open"``, chooses the dialog type
    :param str caption: caption of the dialog
    :param str directory: initial directory
    :param FileTypes file_types: file extensions
    :return: the file name with extension, the description
             of the used filetype and the last path
    :rtype: tuple(str, str, str)
    :raise ValueError: if an invalid value was supplied to ``kind``
    """

    if file_types is None:
        file_types = FileTypes()

    if kind == "save":
        dialog = QFileDialog.getSaveFileName
    elif kind == "open":
        dialog = QFileDialog.getOpenFileName
    else:
        raise ValueError("Unknown value for 'kind': " + kind)

    filename, used_filter = dialog(
        parent,
        caption,
        directory,
        file_types.get_filter())
    new_path = os.path.split(filename)[0]

    if not filename:
        return None, None, new_path

    ext, desc = file_types.get_type(used_filter)
    if not filename.endswith(ext):
        filename += "." + ext
    return filename, desc, new_path
