import os
from typing import Dict, List, Optional, Tuple, Union

from PyQt5.QtWidgets import QFileDialog, QTableWidgetItem, QWidget

_TypeList = Optional[
                Union[
                    List[Tuple[str, str]],
                    Dict[str, str],
                    List[str]
                ]
            ]


class SortableTableWidgetItem(QTableWidgetItem):
    """
    A :class:`QTableWidgetItem` which supports numerical sorting.

    .. automethod:: __init__
    .. automethod:: __lt__
    """

    def __init__(self,
                 text: str) -> None:
        """
        Create a new sortable table widget item.

        :param str text: contents of the item
        :return: nothing
        :rtype: None
        """

        super().__init__(text)

    def __lt__(self,
               other: "SortableTableWidgetItem") -> bool:
        """
        Compare two items.

        :param SortableTableWidgetItem other: item to which self is compared
        :return: True if self is less than other
        :rtype: bool
        """

        key1 = self.text()
        key2 = other.text()
        try:
            return float(key1) < float(key2)
        except ValueError:
            return key1 < key2


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
        "png": ("png", "Portable network graphics"),
        "svg": ("svg", "Scalable vector graphics"),
        "gv": ("gv", "GraphViz DOT"),
        "gexf": ("gexf", "Graph exchange XML format"),
        "": ("", "all files")
    }

    def __init__(self,
                 type_list: _TypeList=None) -> None:
        """
        Create a new :class:`FileTypes` object.

        :param type_list: list of (extension, description) tuples
                          or {ext: description} dict
                          or list of extensions, which are then filtered
                          from ``_default_file_types``
        :type type_list: list(tuple(str, str)) or dict or list(str)
        :return: nothing
        :rtype: None
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
                self.types = []
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

    def get_filter(self) -> str:
        """
        Returns  a file filter suitable for the ``filter`` parameter of
        :class:`~PyQt5.QtWidgets.QFileDialog`.

        :return: a file filter
        :rtype: str
        """

        return ";;".join(self.filters)

    def get_type(self,
                 filefilter: str) -> Tuple[str, str]:
        """
        Returns the extension and description associated with a file filter.

        :param str filefilter: filefilter as returned by the dialog
        :return: the extension and description of the file type
        :rtype: tuple(str, str)
        """

        return self.types[self.filters.index(filefilter)]


def get_filename(parent: QWidget,
                 kind: str="save",
                 caption: str="",
                 directory: str="",
                 file_types: Optional[FileTypes]=None) -> Tuple[Optional[str],
                                                                str]:
    """
    Get a filename by a :class:`QFileDialog`
    and automatically add extensions.

    :param QWidget parent: parent of the dialog
    :param str kind: ``"save"`` or ``"open"``, chooses the dialog type
    :param str caption: caption of the dialog
    :param str directory: initial directory
    :param FileTypes file_types: file extensions
    :return: the file name with extension and the last path
    :rtype: tuple(str, str)
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
        return None, new_path

    ext, _ = file_types.get_type(used_filter)
    if not filename.endswith(ext):
        filename += "." + ext
    return filename, new_path
