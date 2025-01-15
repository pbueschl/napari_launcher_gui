import tifffile
import h5py
import numpy as np
import os
import argparse
import nibabel as nib
import xmltodict
import xml.etree.ElementTree as ET


def read_ome_tiff_image_and_metadata(path_to_ome_tiff_file):
    """
    Opens an OME TIFF file and returns a numpy array of the image data together with a metadata containing dictionary
    :param path_to_ome_tiff_file: path to an OME TIFF file (string)
    :return: image data (numpy.ndarray), metadata (dict)
    """
    with tifffile.TiffFile(path_to_ome_tiff_file) as f:
        # read image data array
        image_data_array = f.asarray()
        # read metadata
        metadata = f.imagej_metadata
        # if channel_names are present in the metadata convert them to a list
        if 'channel_names' in metadata:
            metadata['channel_names'] = eval(metadata['channel_names'])

    # return the metadata
    return image_data_array, metadata


def read_metadata_from_ims_file(path_to_ims_file):
    """
    Reads image data together with relevant metadata from an Imaris ims file and returns the image data together
    with the extracted metadata.

    :param path_to_ims_file: path to the ims file from which the data should be extracted (string)
    :return: image_data_array (numpy.ndarray), metadata_dict (dict)
    """
    # print status message
    print(f'Read "{path_to_ims_file}" ...')
    # open ims file
    with h5py.File(path_to_ims_file, 'r') as f:
        # generate list of available channels
        channel_list = [ch_id for ch_id in f['DataSetInfo'] if
                        ch_id.startswith('Channel') and 'Name' in f['DataSetInfo'][ch_id].attrs.keys()]

        channel_names = []

        # iterate channels
        for channel in channel_list:
            # read channel names
            channel_names.append(f['DataSetInfo'][channel].attrs['Name'].tobytes().decode('ascii', 'ignore'))

        # initialize empty dict for voxel and image sizes
        voxel_size = {}
        image_size = {}
        # read min and max metric coordinates as well as pixel size in all three dimensions and calculate voxel size
        for i, d in [(0, 'X'), (1, 'Y'), (2, 'Z')]:
            # read the highest metrical coordinate
            max_coord = float(f['DataSetInfo']['Image'].attrs[f'ExtMax{i}'].tobytes().decode('ascii', 'decode'))
            # read the lowest metrical coordinate
            min_coord = float(f['DataSetInfo']['Image'].attrs[f'ExtMin{i}'].tobytes().decode('ascii', 'decode'))
            # read the pixel size
            pixel_size = int(f['DataSetInfo']['Image'].attrs[d].tobytes().decode('ascii', 'decode'))

            # calculate metrical voxel size
            voxel_size[d] = (max_coord - min_coord) / pixel_size

            # add dimension size to image size dict
            image_size[d] = pixel_size

    # generate meta data dict
    metadata_dict = {'axes': 'ZCYX',
                     'axes_info': 'ZCYX',
                     'channels': len(channel_names),
                     'slices': image_size['Z'],
                     'hyperstack': True,
                     'mode': 'grayscale',
                     'Channel': [{'Name': s.lower()} for s in channel_names],
                     'image_size': image_size,
                     'voxel_size': voxel_size,
                     'PhysicalSizeX': voxel_size['X'],
                     'PhysicalSizeXUnit': 'µm',
                     'PhysicalSizeY': voxel_size['Y'],
                     'PhysicalSizeYUnit': 'µm',
                     'PhysicalSizeZ': voxel_size['Z'],
                     'PhysicalSizeZUnit': 'µm',
                     'original_file': os.path.basename(path_to_ims_file)}

    # return dictionary with relevant metadata
    return metadata_dict


def read_image_from_ims_file(path_to_ims_file, selected_channel_list=None):
    """
    Reads image data together with relevant metadata from an Imaris ims file and returns the image data together
    with the extracted metadata.

    :param path_to_ims_file: path to the ims file from which the data should be extracted (string)
    :return: image_data_array (numpy.ndarray), metadata_dict (dict)
    """
    # print status message
    print(f'Read "{path_to_ims_file}" ...')
    # open ims file
    with h5py.File(path_to_ims_file, 'r') as f:
        # generate list of available channels
        channel_list = [ch_id for ch_id in f['DataSetInfo'] if
                        ch_id.startswith('Channel') and 'Name' in f['DataSetInfo'][ch_id].attrs.keys()]

        # initialize empty lists for storing the data arrays and names of each channel
        image_data_array = []
        channel_names = []

        if selected_channel_list is None:
            # iterate channels
            for channel in channel_list:
                # read data array
                image_data_array.append(np.array(f['DataSet']['ResolutionLevel 0']['TimePoint 0'][channel]['Data']))
                # read channel names
                channel_names.append(f['DataSetInfo'][channel].attrs['Name'].tobytes().decode('ascii', 'ignore'))
        else:
            # iterate channels
            for channel in channel_list:
                ch_name = f['DataSetInfo'][channel].attrs['Name'].tobytes().decode('ascii', 'ignore')
                if ch_name.lower() in selected_channel_list:
                    # read data array
                    image_data_array.append(np.array(f['DataSet']['ResolutionLevel 0']['TimePoint 0'][channel]['Data']))
                    # read channel names
                    channel_names.append(ch_name)

        # combine 3D arrays into a 4D array
        image_data_array = np.stack(image_data_array)

        # initialize empty dict for voxel and image sizes
        voxel_size = {}
        image_size = {}
        # read min and max metric coordinates as well as pixel size in all three dimensions and calculate voxel size
        for i, d in [(0, 'X'), (1, 'Y'), (2, 'Z')]:
            # read the highest metrical coordinate
            max_coord = float(f['DataSetInfo']['Image'].attrs[f'ExtMax{i}'].tobytes().decode('ascii', 'decode'))
            # read the lowest metrical coordinate
            min_coord = float(f['DataSetInfo']['Image'].attrs[f'ExtMin{i}'].tobytes().decode('ascii', 'decode'))
            # read the pixel size
            pixel_size = int(f['DataSetInfo']['Image'].attrs[d].tobytes().decode('ascii', 'decode'))

            # calculate metrical voxel size
            voxel_size[d] = (max_coord - min_coord) / pixel_size

            # add dimension size to image size dict
            image_size[d] = pixel_size

    # crop array size to image size
    image_data_array = image_data_array[:, :image_size['Z'], :image_size['Y'], :image_size['X']]

    # generate meta data dict
    metadata_dict = {'axes': 'CZYX',
                     'axes_info': 'CZYX',
                     'channels': len(channel_names),
                     'slices': image_size['Z'],
                     'hyperstack': True,
                     'mode': 'grayscale',
                     # 'Channel': [{'ID': f'Channel:0:{ch_id}', 'Name': ch_name.lower()} for ch_id, ch_name in enumerate(channel_names)],
                     'Channel': [{'Name': s.lower()} for s in channel_names],
                     'image_size': image_size,
                     'voxel_size': voxel_size,
                     'PhysicalSizeX': voxel_size['X'],
                     'PhysicalSizeXUnit': 'µm',
                     'PhysicalSizeY': voxel_size['Y'],
                     'PhysicalSizeYUnit': 'µm',
                     'PhysicalSizeZ': voxel_size['Z'],
                     'PhysicalSizeZUnit': 'µm',
                     'original_file': os.path.basename(path_to_ims_file)}

    # return image data array and dictionary with relevant metadata
    return image_data_array, metadata_dict


def get_list_of_channel_names(metadata):
    if type(metadata) is dict:
        if type(metadata['Channel'][0]) is dict:
            list_of_channels = [ch_dict['Name'] for ch_dict in metadata['Channel']]
        else:
            list_of_channels = metadata['Channel']['Name']
    else:
        metadata = xmltodict.parse(metadata)
        list_of_channels = [ch_dict['@Name'] for ch_dict in metadata['OME']['Image']['Pixels']['Channel']]
    return list_of_channels


def subsample_channels(data_array, metadata, desired_channels):
    if type(metadata) is dict:
        if type(metadata['Channel'][0]) is dict:
            dict_of_channels = {ch_dict['Name']: i for i, ch_dict in enumerate(metadata['Channel'])}
        else:
            dict_of_channels = {name: i for i, name in enumerate(metadata['Channel']['Name'])}
        dict_of_new_channels = {dict_of_channels[name]: name for name in desired_channels}
        id_new_ch = np.array(list(dict_of_new_channels.keys()))
        metadata['channels'] = len(id_new_ch)
        metadata['Channel'] = [{'Name': dict_of_new_channels[ch_id]} for ch_id in id_new_ch]

    else:
        metadata = xmltodict.parse(metadata)
        list_of_channels = [ch_dict['@Name'] for ch_dict in metadata['OME']['Image']['Pixels']['Channel']]
        dict_of_channels = {name: i for i, name in enumerate(list_of_channels)}
        dict_of_new_channels = {dict_of_channels[name]: name for name in desired_channels}
        id_new_ch = np.array(list(dict_of_new_channels.keys()))
        metadata['OME']['Image']['Pixels']['@SizeC'] = str(len(id_new_ch))
        metadata['OME']['Image']['Pixels']['Channel'] = [
            {'@ID': f'Channel:0:{i}', '@SamplesPerPixel': '1', '@Name': dict_of_new_channels[ch_id], 'LightPath': None}
            for i, ch_id in enumerate(id_new_ch)]

    # subsample data array
    data_array = data_array[:, id_new_ch, :, :]

    # return data_array and metadata
    return data_array, metadata


def save_nifti_file(image_data_array,
                    metadata_dict,
                    path_to_new_file,
                    default_prefix='image_',
                    default_suffix=('nii', 'nii.gz')):
    """
    Saves passed image and metadata as OME TIFF file. If a directory is passed es destination for saving the image,
    images are stored in consecutive order (dafault_prefix_{number}.default_suffix), if filename is passed image is
    stored utilizing this file name.

    :param image_data_array: array with the image pixel data (numpy.ndarray)
    :param metadata_dict: dictionary with all necessary metadata for the passed image (dict)
    :param path_to_new_ome_file: path to a directory or file name (string)
    :param default_prefix: (string) if not passed default value 'image' is used
    :param default_suffix: (string) if not passed default value 'ome.tif' is used
    """

    # check if path for storing the image ends with the correct file extension (.ome.tif)
    if not path_to_new_file.endswith(default_suffix):
        # if not, consider it as directory and make it if it does not exist
        if not os.path.exists(path_to_new_file):
            os.makedirs(path_to_new_file)
            # print status message
            print(f'Created new directory ({path_to_new_file})!')

    # check if passed path is a directory
    if os.path.isdir(path_to_new_file):
        # get a list of all files in the directory that match the prefix and suffix
        files = [f for f in os.listdir(path_to_new_file) if
                 f.startswith(default_prefix) and f.endswith(default_suffix)]

        # Extract the index numbers from the filenames and sort them
        indices = sorted([int(f[len(default_prefix):-len(default_suffix)]) for f in files])

        # Check for the next available index
        if not indices:
            next_index = 1
        else:
            next_index = indices[-1] + 1
        # Construct the filename for the new file
        new_filename = f"{default_prefix}_{next_index}.{default_suffix}"
        path_to_new_file = os.path.join(path_to_new_file, new_filename)

    # reorder data to (X, Y, Z, C)
    image_data_array = np.transpose(image_data_array, (3, 2, 0, 1))

    # Create an affine matrix
    affine = np.array([[metadata_dict['voxel_size']['X'], 0, 0, 0],
                       [0, metadata_dict['voxel_size']['Y'], 0, 0],
                       [0, 0, metadata_dict['voxel_size']['Z'], 0],
                       [0, 0, 0, 1]])

    # create a NIfTI image
    nifti_img = nib.Nifti1Image(image_data_array, affine)

    # save the NIfTI image
    nifti_img.to_filename(path_to_new_file)

    # print status message
    print(f'Saved file at {path_to_new_file}!')


def save_ome_tiff_file(image_data_array,
                       metadata_dict,
                       path_to_new_ome_file,
                       default_prefix='image_',
                       default_suffix='ome.tif'):
    """
    Saves passed image and metadata as OME TIFF file. If a directory is passed es destination for saving the image,
    images are stored in consecutive order (dafault_prefix_{number}.default_suffix), if filename is passed image is
    stored utilizing this file name.

    :param image_data_array: array with the image pixel data (numpy.ndarray)
    :param metadata_dict: dictionary with all necessary metadata for the passed image (dict)
    :param path_to_new_ome_file: path to a directory or file name (string)
    :param default_prefix: (string) if not passed default value 'image' is used
    :param default_suffix: (string) if not passed default value 'ome.tif' is used
    """

    # check if axis is defined in metadata_dict
    if not 'axes' in metadata_dict:
        # add default axis order to metadata dict
        metadata_dict['axes'] = 'CZYX'

    # check if path for storing the image ends with the correct file extension (.ome.tif)
    if not path_to_new_ome_file.endswith(default_suffix):
        # if not, consider it as directory and make it if it does not exist
        if not os.path.exists(path_to_new_ome_file):
            os.makedirs(path_to_new_ome_file)
            # print status message
            print(f'Created new directory ({path_to_new_ome_file})!')

    # check if passed path is a directory
    if os.path.isdir(path_to_new_ome_file):
        # get a list of all files in the directory that match the prefix and suffix
        files = [f for f in os.listdir(path_to_new_ome_file) if
                 f.startswith(default_prefix) and f.endswith(default_suffix)]

        # Extract the index numbers from the filenames and sort them
        indices = sorted([int(f[len(default_prefix):-len(default_suffix)]) for f in files])

        # Check for the next available index
        if not indices:
            next_index = 1
        else:
            next_index = indices[-1] + 1
        # Construct the filename for the new file
        new_filename = f"{default_prefix}_{next_index}.{default_suffix}"
        path_to_new_ome_file = os.path.join(path_to_new_ome_file, new_filename)

    # save the image data array as ome-tiff file
    with tifffile.TiffWriter(path_to_new_ome_file, bigtiff=True) as tif_writer:
        tif_writer.write(image_data_array,
                         photometric='minisblack',
                         metadata=metadata_dict)
    # print status message
    print(f'Saved file at {path_to_new_ome_file}!')


def convert_ome_xml_metadata_to_dict(xml_metadata):
    # Parse XML
    root = ET.fromstring(xml_metadata)
    pixels = root.find('.//{http://www.openmicroscopy.org/Schemas/OME/2016-06}Pixels')

    # Extracting relevant information
    channel_elements = pixels.findall('.//{http://www.openmicroscopy.org/Schemas/OME/2016-06}Channel')
    channel_names = [channel.get('Name') for channel in channel_elements]

    # Constructing the desired dictionary
    metadata_dict = {
        'axes': 'ZCYX',
        'axes_info': 'ZCYX',
        'channels': len(channel_names),
        'slices': int(pixels.get('SizeZ')),
        'hyperstack': True,
        'mode': 'grayscale',
        'Channel': [{'Name': s.lower()} for s in channel_names],
        'image_size': {'X': int(pixels.get('SizeX')), 'Y': int(pixels.get('SizeY')), 'Z': int(pixels.get('SizeZ'))},
        'voxel_size': {'X': float(pixels.get('PhysicalSizeX')), 'Y': float(pixels.get('PhysicalSizeY')),
                       'Z': float(pixels.get('PhysicalSizeZ'))},
        'PhysicalSizeX': float(pixels.get('PhysicalSizeX')),
        'PhysicalSizeXUnit': pixels.get('PhysicalSizeXUnit'),
        'PhysicalSizeY': float(pixels.get('PhysicalSizeY')),
        'PhysicalSizeYUnit': pixels.get('PhysicalSizeYUnit'),
        'PhysicalSizeZ': float(pixels.get('PhysicalSizeZ')),
        'PhysicalSizeZUnit': pixels.get('PhysicalSizeZUnit'),
    }

    # Return the dictionary
    return metadata_dict
