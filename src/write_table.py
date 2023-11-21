"""Module to write MTG tables"""
import h5py
import numpy as np
import pandas as pd
import datetime
import random
import string

def add_table_header(hdf_file):
    """Common header for MTG tables"""
    hdf_file.attrs["date"] = str(datetime.date.today())
    hdf_file.attrs["target_solver"] = "AVBP"
    hdf_file.attrs["generated_with"] = "MTG"
    hdf_file.attrs["MTG_version"] = 'EGR0D' #__version__
    hdf_file.attrs["HDF5_Version"] = h5py.version.hdf5_version
    hdf_file.attrs["h5py_version"] = h5py.version.version

    return hdf_file

def write_table_hdf(table, chemical_scheme="None", author="None"):
    """Function to write a HDF table for FPI-TTC, NOMAGT, TFLES tools
    with flames data

    Parameters
    ----------
    table : class `GenerateTable`
            Object of class GenerateTable
    chemical_scheme : str, default: `'None'`
            Name of the chemical scheme used
    author : str, default: `'None'`
            Author's name that generated the data
    """
    filename = table.output_folder + "/FLAMES-" + table.mechanism + ".h5"

    phi_list = table.list_sorted_phi

    temperature_list = table.list_sorted_temperature
    pressure_list = table.list_sorted_pressure

    hdf_file = h5py.File(filename, "w")
    hdf_file = add_table_header(hdf_file)
    hdf_file.attrs["file_name"] = filename
    hdf_file.attrs["chemical_scheme"] = chemical_scheme
    hdf_file.attrs["author"] = author

    group = hdf_file.create_group("Data")

    # TODO
    # variables_list = [phi_list, temperature_list, pressure_list]
    # if table.dilution:
    #     variables_list.append(dilution_list)

    # for k, var_list in zip(
    #         list(np.arange(len(variables_list))),
    #         variables_list
    #         ):
    #     for attrs_name in  ['NAME', 'DESCRIPTION', 'NPTS']:
    #         group.attrs[f'DIM{k}_{attrs_name}'] =

    #     group.attrs[f'DIM{k}_MIN'] = var_list[0]
    #     group.attrs[f'DIM{k}_MAX'] = var_list[-1]
    #     group.attrs[f'DIM{k}_VALUES'] = var_list

    if table.phivitiated:
        group.attrs["DIM1_NAME"] = "OXIDIZER_VITIATION"
        group.attrs["DIM1_DESCRIPTION"] = "oxidizer phivitiated"
    else:
        group.attrs["DIM1_NAME"] = "EQUIVALENCE_RATIO"
        group.attrs["DIM1_DESCRIPTION"] = "equivalence ratio"

    group.attrs["DIM1_NPTS"] = len(phi_list)
    group.attrs["DIM1_MIN"] = phi_list[0]
    group.attrs["DIM1_MAX"] = phi_list[-1]
    group.attrs["DIM1_VALUES"] = phi_list

    group.attrs["DIM2_NAME"] = "TEMPERATURE"
    group.attrs["DIM2_DESCRIPTION"] = "temperature of fresh gases"
    group.attrs["DIM2_NPTS"] = len(temperature_list)
    group.attrs["DIM2_VALUES"] = temperature_list
    group.attrs["DIM2_MIN"] = temperature_list[0]
    group.attrs["DIM2_MAX"] = temperature_list[-1]

    group.attrs["DIM3_NAME"] = "PRESSURE"
    group.attrs["DIM3_DESCRIPTION"] = "pressure"
    group.attrs["DIM3_NPTS"] = len(pressure_list)
    group.attrs["DIM3_VALUES"] = pressure_list
    group.attrs["DIM3_MIN"] = pressure_list[0]
    group.attrs["DIM3_MAX"] = pressure_list[-1]

    if table.phivitiated:
        phivitiated_list = table.list_sorted_phivitiated
        group.attrs["DIM4_NAME"] = "EQUIV_RATIO_VITIATED"
        group.attrs["DIM4_DESCRIPTION"] = "equivalence ratio vitiated"
        group.attrs["DIM4_NPTS"] = len(phivitiated_list)
        group.attrs["DIM4_VALUES"] = phivitiated_list
        group.attrs["DIM4_MIN"] = phivitiated_list[0]
        group.attrs["DIM4_MAX"] = phivitiated_list[-1]

    if table.phivitiated:
        phi_ds = group.create_dataset("OXIDIZER_VITIATION", data=phi_list)
    else:
        phi_ds = group.create_dataset("EQUIVALENCE_RATIO", data=phi_list)
    temperature_ds = group.create_dataset("TEMPERATURE",
                                          data=temperature_list)
    pressure_ds = group.create_dataset("PRESSURE", data=pressure_list)

    if table.phivitiated:
        phivitiated_ds = group.create_dataset("EQUIV_RATIO_VITIATED", data=phivitiated_list)

    # Depending on h5py version
    h5py_older_version = False
    try:  # in case h5py>=2.10
        if table.phivitiated:
            phi_ds.make_scale("OXIDIZER_VITIATION")
        else:
            phi_ds.make_scale("EQUIVALENCE_RATIO")
            temperature_ds.make_scale("TEMPERATURE")
            pressure_ds.make_scale("PRESSURE")
            if table.vitiation:
                phivitiated_ds.make_scale("EQUIV_RATIO_VITIATED")
    except BaseException:
        # in case h5py<=2.9
        h5py_older_version = True

    for data_name in table.variable_names:
        # Adding dataset SL_0, TEMPERATURE_BURNT, DELTA_L_0,
        # VOLUMIC_HEAT_RELEASE_MAX, W_YC_MAX, W_FUEL_MAX, Y_*_BURNT
        data_ds = group.create_dataset(data_name, data=table.data[data_name])
        if table.phivitiated:
            data_ds.dims[0].label = "OXIDIZER_VITIATION"
        else:
            data_ds.dims[0].label = "EQUIVALENCE_RATIO"
        data_ds.dims[1].label = "TEMPERATURE"
        data_ds.dims[2].label = "PRESSURE"
        if table.phivitiated:
            data_ds.dims[3].label = "EQUIV_RATIO_VITIATED"
        if h5py_older_version:
            if table.phivitiated:
                data_ds.dims.create_scale(phi_ds, "OXIDIZER_VITIATION")
            else:
                data_ds.dims.create_scale(phi_ds, "EQUIVALENCE_RATIO")
                data_ds.dims.create_scale(temperature_ds, "TEMPERATURE")
                data_ds.dims.create_scale(pressure_ds, "PRESSURE")
            if table.phivitiated:
                data_ds.dims.create_scale(phivitiated_ds, "EQUIV_RATIO_VITIATED")
        data_ds.dims[0].attach_scale(phi_ds)
        data_ds.dims[1].attach_scale(temperature_ds)
        data_ds.dims[2].attach_scale(pressure_ds)
        if table.phivitiated:
            data_ds.dims[3].attach_scale(phivitiated_ds)

    group = hdf_file.create_group("TableID")

    if table.phivitiated:
        group.create_dataset("TableID", data=get_table_identifier(
            phi_list, temperature_list, pressure_list, table.tool, phivitiated_list))
    else:
        group.create_dataset("TableID", data=get_table_identifier(
            phi_list, temperature_list, pressure_list, table.tool))

    group = hdf_file.create_group("Parameters")
    dset = group.create_dataset("STORAGE", (1,), dtype="S09")
    dset[0] = np.string_("ROW-MAJOR")
    group.create_dataset("MIXTURE", data=chemical_scheme)

    if table.phivitiated:
        group.create_dataset("NDIM", data=4)
    else:
        group.create_dataset("NDIM", data=3)

    dset = group.create_dataset("DIM1_NAME", (1,), dtype="S18")
    if table.phivitiated:
        dset[0] = np.string_("OXIDIZER_VITIATION")
    else:
        dset[0] = np.string_("EQUIVALENCE_RATIO")
    dset = group.create_dataset("DIM1_DESCRIPTION", (1,), dtype="S18")
    if table.phivitiated:
        dset[0] = np.string_("oxidizer_phivitiated")
    else:
        dset[0] = np.string_("equivalence_ratio")
    group.create_dataset("DIM1_NPTS", data=len(phi_list))
    dset = group.create_dataset("DIM1_SPACING", (1,), dtype="S05")
    dset[0] = np.string_("INDEX")
    group.create_dataset("DIM1_INDEX", data=phi_list)
    group.create_dataset("NPTS1", data=len(phi_list))
    group.create_dataset("DIM1_MIN", data=phi_list[0])
    group.create_dataset("DIM1_MAX", data=phi_list[-1])
    group.create_dataset("DIM1_VALUES", data=phi_list)

    dset = group.create_dataset("DIM2_NAME", (1,), dtype="S11")
    dset[0] = np.string_("TEMPERATURE")
    dset = group.create_dataset("DIM2_DESCRIPTION", (1,), dtype="S11")
    dset[0] = np.string_("temperature of fresh gases")
    group.create_dataset("DIM2_NPTS", data=len(temperature_list))
    group.create_dataset("NPTS2", data=len(temperature_list))
    group.create_dataset("DIM2_VALUES", data=temperature_list)
    group.create_dataset("DIM2_MIN", data=temperature_list[0])
    group.create_dataset("DIM2_MAX", data=temperature_list[-1])
    dset = group.create_dataset("DIM2_SPACING", (1,), dtype="S05")
    dset[0] = np.string_("INDEX")
    group.create_dataset("DIM2_INDEX", data=temperature_list)

    dset = group.create_dataset("DIM3_NAME", (1,), dtype="S08")
    dset[0] = np.string_("PRESSURE")
    dset = group.create_dataset("DIM3_DESCRIPTION", (1,), dtype="S08")
    dset[0] = np.string_("pressure")
    group.create_dataset("NPTS3", data=len(pressure_list))
    group.create_dataset("DIM3_VALUES", data=pressure_list)
    group.create_dataset("DIM3_MIN", data=pressure_list[0])
    group.create_dataset("DIM3_MAX", data=pressure_list[-1])
    dset = group.create_dataset("DIM3_SPACING", (1,), dtype="S05")
    dset[0] = np.string_("INDEX")
    group.create_dataset("DIM3_INDEX", data=pressure_list)

    if table.phivitiated:
        dset = group.create_dataset("DIM4_NAME", (1,), dtype="S20")
        dset[0] = np.string_("EQUIV_RATIO_VITIATED")
        dset = group.create_dataset("DIM4_DESCRIPTION", (1,), dtype="S20")
        dset[0] = np.string_("phivitiated")
        group.create_dataset("NPTS4", data=len(phivitiated_list))
        group.create_dataset("DIM4_VALUES", data=phivitiated_list)
        group.create_dataset("DIM4_MIN", data=phivitiated_list[0])
        group.create_dataset("DIM4_MAX", data=phivitiated_list[-1])
        dset = group.create_dataset("DIM4_SPACING", (1,), dtype="S05")
        dset[0] = np.string_("INDEX")
        group.create_dataset("DIM4_INDEX", data=phivitiated_list)

    hdf_file.close()
    #logprint("\n\n                 --------------------------------")
    #logprint("  --> WROTE HDF TABLE: " + str(filename))

def get_table_identifier(Phi, T, P, tool: str, D=None) -> str:
    """Function to generate a unique table key ID."""

    def random_identifier(size):
        """Generate the random string of size ``size``"""
        return "".join(
            random.SystemRandom().choice(string.ascii_uppercase + string.digits)
            for _ in range(size)
        )

    # def initials_author_name(author_name):
    #     [first, last] = author_name.split(' ')
    #     return first[0] + last[0]

    # Start of table ID
    table_id = "MTG_" + tool

    # Equivalence ratio range
    if tool in ["FPI-TTC", "NOMAGT", "PAHARC"]:
        first_phi_index = 1
        last_phi_index = -2
    else:
        first_phi_index = 0
        last_phi_index = -1
    table_id += "_Phi" + "{:03{f}}".format(
        Phi[first_phi_index] * 100
        if Phi[first_phi_index] > 1.0
        else int(Phi[first_phi_index] * 100),
        f="g" if Phi[first_phi_index] > 1.0 else "d",
    )
    if len(Phi) > 1:
        table_id += "-" + "{:03{f}}".format(
            Phi[last_phi_index] * 100
            if Phi[last_phi_index] > 1.0
            else int(Phi[last_phi_index] * 100),
            f="g" if Phi[last_phi_index] > 1.0 else "d",
        )

    # Temperature range
    table_id += f"_T{int(T[0]):04d}"
    if len(T) > 1:
        table_id += f"-{int(T[-1]):04d}"

    # Pressure range
    table_id += f"_P{int(P[0]/1000):04d}"
    if len(P) > 1:
        table_id += f"-{int(P[-1]/1000):04d}"

    if D is not None:
        table_id += "_Dil{D[0]:03d}-{D[-1]:03d}"

    # Random key
    table_id += "_id" + random_identifier(5)

    return table_id