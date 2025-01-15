import os
import napari
from PyQt5 import QtWidgets, QtGui, QtCore
from qtpy.QtWidgets import QPushButton, QVBoxLayout, QWidget, QFileDialog
import tifffile as tif
from utils.image_file_io import convert_ome_xml_metadata_to_dict, read_image_from_ims_file
import nibabel as nib


class NapariLauncher(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.files_to_load = []  # Store files with associated settings
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Napari Launcher")
        self.setGeometry(100, 100, 800, 400)

        layout = QtWidgets.QVBoxLayout()

        # File selection
        self.file_label = QtWidgets.QLabel("File to load:")
        self.file_input = QtWidgets.QLineEdit()
        self.file_button = QtWidgets.QPushButton("Browse")
        self.file_button.clicked.connect(self.browse_file)

        file_layout = QtWidgets.QHBoxLayout()
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(self.file_button)

        # Brush size input and data type selection
        self.brush_size_label = QtWidgets.QLabel("Maximum brush size:")
        self.brush_size_input = QtWidgets.QSpinBox()
        self.brush_size_input.setRange(1, 5000)
        self.brush_size_input.setValue(1000)

        self.data_type_label = QtWidgets.QLabel("Data type:")
        self.data_type_combo = QtWidgets.QComboBox()
        self.data_type_combo.addItems(["single_class_label", "multi_class_label", "image"])
        self.data_type_combo.currentIndexChanged.connect(self.toggle_brush_size)

        brush_data_layout = QtWidgets.QHBoxLayout()
        brush_data_layout.addWidget(self.brush_size_label)
        brush_data_layout.addWidget(self.brush_size_input)
        brush_data_layout.addWidget(self.data_type_label)
        brush_data_layout.addWidget(self.data_type_combo)

        # Add file button
        self.add_file_button = QtWidgets.QPushButton("Add File")
        self.add_file_button.clicked.connect(self.add_file)

        # Spacer (fixed height equal to button height)
        spacer = QtWidgets.QSpacerItem(20, self.add_file_button.sizeHint().height(), QtWidgets.QSizePolicy.Minimum,
                                       QtWidgets.QSizePolicy.Fixed)

        # File table
        self.file_table = QtWidgets.QTableWidget()
        self.file_table.setColumnCount(3)
        self.file_table.setHorizontalHeaderLabels(["File Name", "Brush Size", "Data Type"])
        self.file_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.file_table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
        self.file_table.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeToContents)
        self.file_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.file_table.customContextMenuRequested.connect(self.show_context_menu)

        # Launch button
        self.launch_button = QtWidgets.QPushButton("Launch Napari")
        self.launch_button.clicked.connect(self.launch_napari)

        # Adding widgets to layout
        layout.addWidget(self.file_label)
        layout.addLayout(file_layout)
        layout.addLayout(brush_data_layout)
        layout.addWidget(self.add_file_button)
        layout.addSpacerItem(spacer)
        layout.addWidget(self.file_table)
        layout.addWidget(self.launch_button)
        self.setLayout(layout)

    def browse_file(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select File", "",
                                                             "All Supported Files (*.tif *.tiff *.ome.tif *.nii *.nii.gz *.ims)")
        if file_path:
            self.file_input.setText(file_path)

    def toggle_brush_size(self):
        if self.data_type_combo.currentText() == "image":
            self.brush_size_input.setEnabled(False)
        else:
            self.brush_size_input.setEnabled(True)

    def add_file(self):
        file_path = self.file_input.text()
        max_brush_size = self.brush_size_input.value()
        data_type = self.data_type_combo.currentText()

        if not file_path:
            QtWidgets.QMessageBox.warning(self, "Error", "Please select a file to add.")
            return

        file_info = {
            "file_path": file_path,
            "max_brush_size": max_brush_size if data_type != "image" else None,
            "data_type": data_type
        }
        self.files_to_load.append(file_info)

        # Update table
        row_position = self.file_table.rowCount()
        self.file_table.insertRow(row_position)

        # Add file name with tooltip for full path
        file_name_item = QtWidgets.QTableWidgetItem(os.path.basename(file_path))
        file_name_item.setToolTip(file_path)
        self.file_table.setItem(row_position, 0, file_name_item)

        # Add brush size
        self.file_table.setItem(row_position, 1,
                                QtWidgets.QTableWidgetItem(str(max_brush_size) if data_type != "image" else "N/A"))

        # Add data type
        self.file_table.setItem(row_position, 2, QtWidgets.QTableWidgetItem(data_type))

        self.file_input.clear()

    def show_context_menu(self, position):
        menu = QtWidgets.QMenu()

        remove_action = menu.addAction("Remove File")
        edit_action = menu.addAction("Edit File Info")
        edit_scale_action = menu.addAction("Edit Custom Scale")

        action = menu.exec_(self.file_table.viewport().mapToGlobal(position))

        selected_row = self.file_table.currentRow()
        if action == remove_action:
            self.remove_file(selected_row)
        elif action == edit_action:
            self.edit_file(selected_row)
        elif action == edit_scale_action:
            self.edit_custom_scale(selected_row)

    def remove_file(self, row):
        if row >= 0:
            self.file_table.removeRow(row)
            del self.files_to_load[row]

    def edit_file(self, row):
        if row >= 0:
            file_info = self.files_to_load[row]

            new_file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Edit File", file_info["file_path"],
                                                                     "TIFF Files (*.tif *.tiff)")
            if new_file_path:
                file_info["file_path"] = new_file_path

            new_brush_size, ok = QtWidgets.QInputDialog.getInt(self, "Edit Brush Size", "Maximum Brush Size:",
                                                               file_info["max_brush_size"] or 1, 1, 5000)
            if ok:
                file_info["max_brush_size"] = new_brush_size

            new_data_type, ok = QtWidgets.QInputDialog.getItem(self, "Edit Data Type", "Data Type:",
                                                               ["single_class_label", "multi_class_label", "image"],
                                                               editable=False)
            if ok:
                file_info["data_type"] = new_data_type

            # Update table
            self.file_table.item(row, 0).setText(os.path.basename(file_info["file_path"]))
            self.file_table.item(row, 0).setToolTip(file_info["file_path"])
            self.file_table.item(row, 1).setText(
                str(file_info["max_brush_size"]) if file_info["data_type"] != "image" else "N/A")
            self.file_table.item(row, 2).setText(file_info["data_type"])

    def edit_custom_scale(self, row):
        if row >= 0:
            file_info = self.files_to_load[row]

            if "custom_scale" in file_info.keys():
                current_scale = file_info["custom_scale"]
            else:
                current_scale = None
            current_scale_str = ",".join(map(str, current_scale)) if current_scale else ""

            new_custom_scale, ok = QtWidgets.QInputDialog.getText(self, "Edit Custom Scale",
                                                                  f"Enter custom scale (comma-separated Z,Y,X):",
                                                                  text=current_scale_str)
            if ok and new_custom_scale:
                file_info["custom_scale"] = [float(v) for v in new_custom_scale.split(",")]
            else:
                file_info["custom_scale"] = None

    def get_data_array_scale_and_channel_axis(self, file_path):

        # read file extension
        file_extension = '.'.join(os.path.basename(file_path).split('.')[1:])
        # TIFF files
        if file_extension == 'ome.tif' or file_extension == 'tif':
            with tif.TiffFile(file_path) as f:
                data_array = f.asarray()

            # try to read some metadata
            try:
                # read metadata
                metadata = f.ome_metadata
                if metadata:
                    # convert xml string to tif
                    metadata = convert_ome_xml_metadata_to_dict(metadata)
                    voxel_size = metadata['voxel_size']
                else:
                    metadata = f.imagej_metadata
                    voxel_size = eval(metadata['voxel_size'])
                scale = (voxel_size['Z'], voxel_size['Y'], voxel_size['X'])

                # check if it is a multichannel image
                axes_info = metadata['axes_info']
                # check if axes info is of same length as array dimensions
                if len(axes_info) > len(data_array.shape):
                    axes_info = axes_info.replace('C', '')
                    metadata['axes_info'] = axes_info
                # add axes key to metadata
                metadata['axes'] = axes_info

                if 'C' in axes_info:
                    channel_axis = axes_info.index('C')
                else:
                    channel_axis = None

            except:
                # if not possible return default scale and channel axis
                scale = scale = [1 for i in data_array.shape]
                channel_axis = None
        # NIFTI Files
        elif file_extension == 'nii' or file_extension == 'nii.gz':
            f = nib.load(file_path)
            data_array = f.get_fdata()
            affine_matrix = f.affine
            scale = [affine_matrix[i, i] for i in range(len(data_array.shape))][:3]
            channel_axis = None
            if len(data_array.shape) > 3:
                channel_axis = len(data_array.shape) - 1
            metadata = {'axes': 'ZYX',
                     'axes_info': 'ZYX'}

        # Imaris files
        elif file_extension == 'ims':
            # read image and metadata from ims file
            data_array, metadata = read_image_from_ims_file(file_path)
            voxel_size = metadata['voxel_size']
            scale = (voxel_size['Z'], voxel_size['Y'], voxel_size['X'])
            axes_info = metadata['axes_info']
            channel_axis = axes_info.index('C')

        else:
            # throw error that file_type is not supported
            raise ValueError(f"Unsupported file type: {file_extension}")

        return data_array, scale, channel_axis, metadata

    def launch_napari(self):
        if not self.files_to_load:
            QtWidgets.QMessageBox.warning(self, "Error", "No files to load. Please add files first.")
            return

        viewer = napari.Viewer()

        # Create a QWidget for the button
        custom_widget = QWidget()
        layout = QVBoxLayout()

        # Create the save button
        save_button = QPushButton("Save As TIFF")
        save_button.clicked.connect(lambda: self.custom_tif_save_function(viewer))

        # Add the button to the layout and set it in the widget
        layout.addWidget(save_button)
        custom_widget.setLayout(layout)

        # Add the custom widget to the napari GUI
        viewer.window.add_dock_widget(custom_widget, name="Custom Tiff Save As", area="right")

        for file_info in self.files_to_load:
            file_path = file_info["file_path"]
            max_brush_size = file_info["max_brush_size"]
            data_type = file_info["data_type"]

            data_array, scale, channel_axis, metadata = self.get_data_array_scale_and_channel_axis(file_path)
            # check if custom scale was defined
            if "custom_scale" in file_info.keys():
                scale = file_info["custom_scale"]
            axis_info = 'ZYX' if len(scale) == 3 else 'YX'
            if data_type == 'multi_class_label':
                viewer.add_labels(data_array, name=os.path.basename(file_path), scale=scale, metadata=metadata)
                labels_layer = viewer.layers[os.path.basename(file_path)]
                labels_layer.brush_size = max_brush_size
            elif data_type == 'single_class_label':
                data_array = data_array.astype(bool)
                viewer.add_labels(data_array, name=os.path.basename(file_path), scale=scale, metadata=metadata)
                labels_layer = viewer.layers[os.path.basename(file_path)]
                #labels_layer.metadata = {"DimensionOrder": axis_info}
                labels_layer.brush_size = max_brush_size
            elif data_type == 'image':
                # viewer.add_image(data_array, name=os.path.basename(file_path), scale=scale, channel_axis=channel_axis, metadata=metadata)
                viewer.add_image(data_array, name=os.path.basename(file_path), scale=scale, metadata=metadata)

        napari.run()

    def custom_tif_save_function(self ,viewer):
        # Open a file dialog to select save location
        save_path, _ = QFileDialog.getSaveFileName(
            caption="Save Selected Layer As TIFF File",
            filter="All Files (*);;Tiff Files (*.tif)"
        )

        if save_path:
            # Ensure the file name ends with the correct extension
            if not save_path.endswith(".tif") and not save_path.endswith(".ome.tif"):
                # Default to `.tif` if no valid extension is provided
                save_path += ".ome.tif"
            # Perform the saving operation
            print(f"Saving to: {save_path}")

            # save combined mask to tif file
            tif.imwrite(
                save_path,
                viewer.layers.selection.active.data,
                metadata=viewer.layers.selection.active.metadata,
                # description=ome_metadata,  # OME-XML metadata
                bigtiff=True  # Force BigTIFF format
            )


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    launcher = NapariLauncher()
    launcher.show()
    app.exec_()
