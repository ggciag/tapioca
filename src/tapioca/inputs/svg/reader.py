'''
Build a Mandyoc interfaces.txt from an SVG drawing (e.g. made in Inkscape).

Format: SVG 1.1 files saved by Inkscape (tested with Inkscape 1.4, see SVG_VERSION and INKSCAPE_VERSION). The
labels, titles and descriptions of the shapes are the Inkscape object properties (inkscape:label, <title> and
<desc>), and the shapes are read with svgelements, so curves, transforms and groups are supported. SVGScenario
and tapioca-svg-interfaces warn when a file declares another SVG version or was saved by another version of
Inkscape (svg_versions), since a change in the format could change how the drawing is read.

Each filled shape is a body of rock described by its object properties
(Inkscape: Object > Object Properties):

    Label        name of the unit (layer_label). Falls back to the id.
    Title        optional, used in the legend of the plot.
    Description  first line: rheology (a material of tapioca.DEFAULT_MATERIALS:
                 wet olivine, dry olivine, wet quartz, wet anorthite, linear salt,
                 dry quartz, basalt, plagioclase, air); next lines: key = value
                 with the layer properties, e.g.

                     wet quartz
                     C = 1.0
                     rho = 2800.0
                     H = 0.19e-6/2800.0

                 Keys: C, rho, H (required), and A, n, Q, V, k, weakening_seed,
                 cohesion_min, cohesion_max, friction_angle_min and
                 friction_angle_max (optional, they override the values of the
                 rheology in DEFAULT_MATERIALS). The names of
                 tapioca.INTERFACES_PARAMETERS (e.g. density, scale_factor,
                 radiogenic_heat) are also accepted.

                 Optional: h = 120 km (or h = 120e3 m), the reference
                 thickness of the unit, so the drawing does not need to be in
                 vertical scale. See rescale_thickness:
                   - layers (units that reach both sides of the domain): their most
                     common drawn thickness becomes h and the same factor is
                     applied in every column, so thinned regions keep their
                     drawn proportion. Layers without h keep the drawn thickness;
                   - bodies (seeds, basins, weak zones, ...): without h they
                     follow the scaling of the unit they are drawn over; with
                     h they get that thickness and the unit around them
                     absorbs the difference, so the units below do not move;
                   - the bottom unit (e.g. asthenosphere) fills the rest of
                     the domain, so it can not have h.

Shapes with the same label are the same unit, so the description only needs
to be written once. A shape without description inherits the one of the group
(e.g. an Inkscape layer) that contains it. Shapes without any description are
ignored (annotations, arrows, ...). The fill colour is used only in the plot.

Domain: a shape labelled "setup" (usually a rectangle without fill, e.g. in its
own Inkscape layer named "domain") gives the model domain. Its bounding box is mapped to the
domain, top-left corner = (x=0, z=0), bottom-right = (Lx, -Lz), and its
description gives the dimensions and the resolution of the Mandyoc grid:

                     Lx = 1600 km
                     Lz = 400 km
                     dx = 2 km        (or Nx = 801)
                     dz = 2 km        (or Nz = 201)
                     min_thickness = 1 km   (optional, see remove_slivers; 0 km keeps the slivers)

It may also give the settings of the numerical experiment, which are returned
as text in setup["experiment"] and read by the
scenario builder (tapioca.inputs.svg.SVGScenario): regime (rifting, compression, accordion, quiescence), velocity
(e.g. 1 cm/yr), moving_boundary (left, right or symmetric), time_max (e.g. 100 Myr),
quiescence_time (e.g. 0 Myr), magmatism (on or off), magmatism_extraction (on or off)
magmatism_Ts0 (e.g. 1080, in oC) and preset_initial_temperature (folder or .nc file of a
preliminary scenario whose final thermal state is used as initial temperature).

Without a setup shape the SVG page is the domain and Lx, Lz, Nx must be given.
x and z are scaled independently, so the drawing may be vertically
exaggerated. Shapes are painted in document order (later shapes cover earlier
ones), as Inkscape renders them. Points not covered by any shape take the unit
of the nearest covered point (see fill_gaps), and thin layers left by shapes
that overlap a little are removed (see remove_slivers), with warnings.

Since Mandyoc layers are bounded by single-valued interfaces z(x), the script
reads, for every column x, the sequence of units from bottom to top and finds a
single ordered stack containing all of them (shortest common supersequence).
Units absent in a column get zero thickness there. A seed inside the lower
crust becomes [..., lower crust, seed, lower crust, ...]; away from the seed its
two interfaces stay at the depth of the top of the seed, so that Mandyoc, which
interpolates the interfaces between the columns, draws the seed as a box and not
as a funnel up to the top of the lower crust. Each slot of the stack is a
MandyocLayer whose interface is the depth of its base.

Usage
-----
    tapioca-svg-interfaces scenario.svg -o interfaces.txt --plot setup.png
    tapioca-svg-interfaces scenario.svg --dx 1     (override the horizontal resolution of the SVG)

or, from another script:

    from tapioca.inputs.svg import read_setup, layers_from_svg
    setup = read_setup("scenario.svg")
    layers = layers_from_svg("scenario.svg", setup["Lx"], setup["Lz"], setup["Nx"])
    builder = layers_to_builder(layers, setup["Lx"], setup["Lz"], Nz=setup["Nz"])  # tapioca.inputs.MandyocBuilder

Requires: numpy, matplotlib, svgelements (pip install svgelements)
'''

import argparse
import ast
import copy
import operator
import os
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import Patch
import svgelements as se

from ..._variables import DEFAULT_MATERIAL, INTERFACES_PARAMETERS, UNITS
from ..._aux_functions import split_unit
from ..generic_input_classes import MandyocBuilder
from ..checks import SVGInputError, SVGChecks
from .materials import MandyocLayer, get_rheology

__all__ = ["SVGInputError", "SVG_VERSION", "INKSCAPE_VERSION", "svg_versions", "read_setup", "read_units", "layers_from_svg", "layers_to_builder", "write_interfaces",
           "layer_indices", "slot_thickness", "layer_thickness",
           "plot_layers", "points_from_spacing", "EXPERIMENT_KEYS", "SETUP_LABEL"]


SVG_NS = "{http://www.w3.org/2000/svg}"
INKSCAPE_NS = "{http://www.inkscape.org/namespaces/inkscape}"
INKSCAPE_LABEL = INKSCAPE_NS + "label"
SHAPE_TAGS = {"rect", "path", "circle", "ellipse", "polygon", "polyline"}
SVG_VERSION = "1.1"  # version of the SVG format read by tapioca
INKSCAPE_VERSION = "1.4"  # version of Inkscape (major.minor) used to draw and test the SVG scenarios
SETUP_LABEL = "setup"  # label (or id) of the shape that gives the model domain, the resolution and the experiment
MIN_THICKNESS = 1.0e3  # default min_thickness [m]: thinner pieces of units drawn thicker elsewhere are removed (see remove_slivers)
EXPERIMENT_KEYS = ("regime", "velocity", "moving_boundary", "time_max", "quiescence_time",
                   "magmatism", "magmatism_extraction", "magmatism_Ts0", "preset_initial_temperature",
                   "temperature_solver", "water_layer", "magmatic_layer")  # settings of the experiment in the domain description

# key in the description -> parameter of interfaces.txt: the keys of DEFAULT_MATERIAL (rho, C, H, ...) and their
# names in INTERFACES_PARAMETERS (density, scale_factor, radiogenic_heat, ...)
PROPERTY_KEYS = {**{key: key for key in DEFAULT_MATERIAL}, **{name: key for key, name in INTERFACES_PARAMETERS.items()}}

###############################################################################
# Functions
###############################################################################

# Arithmetic operators accepted in the values of the descriptions (e.g. H = 1.67e-6 / 2700.0). eval_number evaluates
# the syntax tree of the value with them, never with eval(), so a description can not run code.
_OPERATORS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
              ast.Div: operator.truediv, ast.Pow: operator.pow, ast.USub: operator.neg, ast.UAdd: operator.pos}

def eval_number(expr):
    '''Evaluate a numeric expression such as "1.67e-6 / 2700.0" (only numbers and + - * / **).'''
    def _eval(node):
        if isinstance(node, ast.Expression):
            return _eval(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.BinOp) and type(node.op) in _OPERATORS:
            return _OPERATORS[type(node.op)](_eval(node.left), _eval(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in _OPERATORS:
            return _OPERATORS[type(node.op)](_eval(node.operand))
        raise ValueError(f"invalid number: {expr!r}")
    return _eval(ast.parse(expr.strip(), mode="eval"))

def parse_quantity(value, quantity, where, name, default_unit=None, positive=False):
    '''
    Read a number with a unit of tapioca.UNITS[quantity] (e.g. "120 km", "1 cm/yr", "300 Myr") and return it
    in the base unit of the quantity (m, yr, cm/yr or oC). The number may be an expression (see eval_number).
    Without unit, default_unit is used; if it is None, the unit is required.
    '''
    units = UNITS[quantity]
    number, factor = split_unit(value, quantity)
    if factor is None:
        if default_unit is None:
            raise SVGInputError(f"{where}: give the unit of {name} ({' or '.join(units)}), e.g. {name} = 1 {next(iter(units))}")
        factor = units[default_unit]
    try:
        result = eval_number(number) * factor
    except (ValueError, SyntaxError):
        raise SVGInputError(f"{where}: could not read the value of {name} = {value}. Units: {', '.join(units)}")
    if positive and result <= 0:
        raise SVGInputError(f"{where}: {name} must be positive.")
    return result


def parse_thickness(value, where, name="h"):
    '''Read a length such as "120 km" or "120e3 m" and return it in m.'''
    return parse_quantity(value, "length", where, name, positive=True)


def points_from_spacing(length, spacing, name):
    '''Number of grid points for a grid spacing; the spacing must divide the length.'''
    intervals = length / spacing
    if abs(intervals - round(intervals)) > 1e-6:
        raise SVGInputError(f"{name} = {spacing / 1e3:g} km does not divide the domain ({length / 1e3:g} km).")
    return int(round(intervals)) + 1

def is_setup_shape(el, label=SETUP_LABEL):
    '''True if el is a shape labelled (or with id) label, by default the setup shape.'''
    return (el.tag.replace(SVG_NS, "") in SHAPE_TAGS
            and label in (el.get(INKSCAPE_LABEL), el.get("id")))

## all this things could be used to create the MandyocBuilder (next implementation: create the abstract class `GeometricDomain`) 
def svg_versions(svg_file):
    '''
    Version of the SVG format (attribute version of <svg>) and version of Inkscape that saved the file
    (inkscape:version, e.g. "1.4.4 (dcaf3e7, 2026-05-05)"), or None when the file does not declare them.
    '''
    root = ET.parse(svg_file).getroot()
    return root.get("version"), root.get(INKSCAPE_NS + "version")


def read_setup(svg_file):
    '''
    Read the shape labelled "setup" (see the module docstring).
    Return dict(id, Lx, Lz, Nx, Nz, experiment and min_thickness if given) with lengths in m, or None if the
    SVG has no setup shape.
    experiment is a dict with the settings of the numerical experiment given in the description (as text).
    '''
    root = ET.parse(svg_file).getroot()
    label = SETUP_LABEL
    shapes = [el for el in root.iter() if is_setup_shape(el, label)]
    if not shapes:
        return None
    if len(shapes) > 1:
        raise SVGInputError(f"The SVG has {len(shapes)} shapes labelled '{label}'; keep only one.")
    el = shapes[0]
    where = f"Setup shape (id={el.get('id')})"
    if el.get("id") is None:
        raise SVGInputError(f"{where}: the setup shape has no id. Save the file with Inkscape to create it.")
    desc = el.find(SVG_NS + "desc")
    if desc is None or not (desc.text or "").strip():
        raise SVGInputError(f"{where}: write Lx, Lz and the resolution (dx and dz, or Nx and Nz) in its description.")

    values, experiment = {}, {}
    for line in desc.text.splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise SVGInputError(f"{where}: unexpected line '{line}' in the description.")
        key, value = (s.strip() for s in line.split("=", 1))
        if key == "min_thickness" and re.fullmatch(r"0*\.?0*\s*(km|m)?", value):
            values[key] = 0.0  # keep the slivers
        elif key in ("Lx", "Lz", "dx", "dz", "min_thickness"):
            values[key] = parse_thickness(value, where, key)
        elif key in ("Nx", "Nz"):
            try:
                values[key] = int(value)
            except ValueError:
                raise SVGInputError(f"{where}: {key} must be an integer.")
        elif key in EXPERIMENT_KEYS:
            experiment[key] = value
        else:
            raise SVGInputError(f"{where}: unknown key '{key}'. Use Lx, Lz, dx, dz (or Nx, Nz), min_thickness and {', '.join(EXPERIMENT_KEYS)}.")

    for key in ("Lx", "Lz"):
        if key not in values:
            raise SVGInputError(f"{where}: missing {key} in the description.")
    setup = dict(id=el.get("id"), Lx=values["Lx"], Lz=values["Lz"], experiment=experiment)
    if "min_thickness" in values:
        setup["min_thickness"] = values["min_thickness"]
    for d, N, L in (("dx", "Nx", "Lx"), ("dz", "Nz", "Lz")):
        if d in values:
            setup[N] = points_from_spacing(values[L], values[d], d)
            if N in values and values[N] != setup[N]:
                raise SVGInputError(f"{where}: {d} and {N} do not agree; give only one of them.")
        elif N in values:
            setup[N] = values[N]
        else:
            raise SVGInputError(f"{where}: give the resolution {d} (or {N}) in the description.")
    return setup

## If MandyocLayer is really redundant, this code have to change
def parse_description(text, where):
    '''
    Parse the description of a shape into
    (rheology (LithologicalUnit), {parameter of DEFAULT_MATERIAL: value}, reference thickness in m or None).
    '''
    rheology, props, thickness = None, {}, None
    for line in text.splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        if "=" not in line:
            if rheology is not None:
                raise SVGInputError(f"{where}: unexpected line '{line}' in the description.")
            rheology = line
            continue
        key, value = (s.strip() for s in line.split("=", 1))
        if key.lower() == "rheology":
            rheology = value
            continue
        if key == "h":
            thickness = parse_thickness(value, where)
            continue
        if key not in PROPERTY_KEYS:
            raise SVGInputError(f"{where}: unknown property '{key}'. Use {', '.join(DEFAULT_MATERIAL)}.")
        try:
            props[PROPERTY_KEYS[key]] = eval_number(value)
        except (ValueError, SyntaxError):
            raise SVGInputError(f"{where}: could not read the value of {key} = {value}")

    if rheology is None:
        raise SVGInputError(f"{where}: the first line of the description must be the rheology (e.g. wet quartz).")
    try:
        rheology = get_rheology(rheology)
    except ValueError as err:
        raise SVGInputError(f"{where}: {err}")
    missing = [s for s in ["C", "rho", "H"] if s not in props]
    if missing:
        raise SVGInputError(f"{where}: missing {', '.join(missing)} in the description.")
    return rheology, props, thickness


def read_units(svg_file):
    '''
    Read label, title and description of the shapes.
    Return {shape id: unit label} and {unit label: dict(rheology, props, title)}.
    '''
    root = ET.parse(svg_file).getroot()
    parent = {child: el for el in root.iter() for child in el}

    def child_text(el, tag):
        node = el.find(SVG_NS + tag)
        return node.text.strip() if node is not None and node.text and node.text.strip() else None

    def label_of(el):
        return el.get(INKSCAPE_LABEL) or el.get("id")

    shape_unit, units = {}, {}
    pending = []  # shapes with a label but no description of their own

    for el in root.iter():
        if el.tag.replace(SVG_NS, "") not in SHAPE_TAGS or is_setup_shape(el):
            continue
        # nearest element (itself or a group above it) with a description
        node = el
        while node is not None and child_text(node, "desc") is None:
            node = parent.get(node)

        if node is None:
            if el.get(INKSCAPE_LABEL):
                pending.append(el)
            continue

        label = label_of(node)
        where = f"Shape '{label}' (id={node.get('id')})"
        rheology, props, thickness = parse_description(child_text(node, "desc"), where)
        unit = dict(rheology=rheology, props=props, thickness=thickness, title=child_text(node, "title"))
        if label in units and (units[label]["rheology"], units[label]["props"], units[label]["thickness"]) != (rheology, props, thickness):
            raise SVGInputError(f"{where}: another shape with the same label has a different description.")
        units.setdefault(label, unit)
        shape_unit[el.get("id")] = label

    for el in pending:
        label = el.get(INKSCAPE_LABEL)
        if label in units:
            shape_unit[el.get("id")] = label
        else:
            print(f"Warning: shape '{label}' (id={el.get('id')}) has no description and was ignored.")

    if None in shape_unit:
        raise SVGInputError("Some shapes with description have no id. Save the file with Inkscape to create them.")
    return shape_unit, units


def shape_to_mplpath(shape, n_curve=32):
    '''Convert an svgelements shape (already transformed) into a matplotlib Path.'''
    path = se.Path(shape)
    verts, codes = [], []
    for seg in path:
        if isinstance(seg, se.Move):
            verts.append((seg.end.x, seg.end.y))
            codes.append(MplPath.MOVETO)
        elif isinstance(seg, se.Close):
            verts.append((seg.end.x, seg.end.y))
            codes.append(MplPath.CLOSEPOLY)
        elif isinstance(seg, se.Line):
            verts.append((seg.end.x, seg.end.y))
            codes.append(MplPath.LINETO)
        else:  # curves and arcs: sample them
            for p in seg.npoint(np.linspace(0, 1, n_curve)[1:]):
                verts.append((p[0], p[1]))
                codes.append(MplPath.LINETO)
    if not verts:
        return None
    return MplPath(np.array(verts), codes)


def rasterize(svg_file, shape_unit, unit_names, Lx, Lz, Nx, dz, setup_id=None):
    '''
    Sample the SVG on a (Nzs, Nx) grid and return the unit index of each point (-1 = empty)
    and the unit it covers (the one that would be there without it, -1 = none).
    Also return the fill colour of each unit (from its first shape).
    '''
    svg = se.SVG.parse(svg_file, reify=True)

    # SVG region mapped to the domain: page, or the bounding box of an element
    if setup_id:
        ref = svg.get_element_by_id(setup_id)
        if ref is None:
            raise SVGInputError(f"Element with id '{setup_id}' not found in the SVG.")
        x0, y0, x1, y1 = ref.bbox(with_stroke=False)
    else:
        x0, y0, x1, y1 = 0.0, 0.0, float(svg.width), float(svg.height)

    Nzs = int(round(Lz / dz))
    dz = Lz / Nzs
    x = np.linspace(0, Lx, Nx)
    zdepth = (np.arange(Nzs) + 0.5) * dz  # depth of samples (positive downwards)

    # sample points in SVG coordinates
    XS, YS = np.meshgrid(x0 + x / Lx * (x1 - x0), y0 + zdepth / Lz * (y1 - y0))
    pts = np.column_stack([XS.ravel(), YS.ravel()])

    grid = -np.ones(XS.shape, dtype=int)
    under = -np.ones(XS.shape, dtype=int)  # unit covered by the one painted at each point (-1 = none)
    colors = {}
    for el in svg.elements():
        if not isinstance(el, se.Shape) or el.id not in shape_unit:
            continue
        if el.values.get("display") == "none" or el.values.get("visibility") == "hidden":
            continue
        iunit = unit_names.index(shape_unit[el.id])
        if el.fill is not None and el.fill.value is not None:
            colors.setdefault(iunit, el.fill.hex[:7])
        mpath = shape_to_mplpath(el)
        if mpath is None:
            continue
        bx0, by0 = mpath.vertices.min(axis=0)
        bx1, by1 = mpath.vertices.max(axis=0)
        sel = (pts[:, 0] >= bx0) & (pts[:, 0] <= bx1) & (pts[:, 1] >= by0) & (pts[:, 1] <= by1)
        if not sel.any():
            continue
        inside = np.zeros(len(pts), dtype=bool)
        inside[sel] = mpath.contains_points(pts[sel])
        covered = inside & (grid.ravel() != iunit)
        under.ravel()[covered] = grid.ravel()[covered]
        grid.ravel()[inside] = iunit

    return x, zdepth, dz, grid, under, colors

def scs(a, b):
    '''
    Shortest common supersequence of two sequences: the shortest sequence that contains both, in order.
    build_stack uses it to find one stack of layers (bottom to top) that fits every column of the drawing:
    the columns [asthenosphere, mantle, crust, air] and [asthenosphere, mantle, crust, seed, crust, air] give
    [asthenosphere, mantle, crust, seed, crust, air], and the first column gets zero thickness for the seed and
    the second crust slot. Mandyoc needs this single stack because each interface is one z(x) over the domain.
    '''
    m, n = len(a), len(b)
    L = np.zeros((m + 1, n + 1), dtype=int)
    L[m, :] = np.arange(n, -1, -1)
    L[:, n] = np.arange(m, -1, -1)
    for i in range(m - 1, -1, -1):
        for j in range(n - 1, -1, -1):
            L[i, j] = 1 + (L[i + 1, j + 1] if a[i] == b[j] else min(L[i + 1, j], L[i, j + 1]))
    out, i, j = [], 0, 0
    while i < m and j < n:
        if a[i] == b[j]:
            out.append(a[i]); i += 1; j += 1
        elif L[i + 1, j] <= L[i, j + 1]:
            out.append(a[i]); i += 1
        else:
            out.append(b[j]); j += 1
    return out + list(a[i:]) + list(b[j:])


def embed(seq, stack):
    '''Slot index in stack of each element of seq (leftmost subsequence match).'''
    slots, j = [], 0
    for m in seq:
        while stack[j] != m:
            j += 1
        slots.append(j)
        j += 1
    return slots


def nearest_covered(covered):
    '''Index of the nearest True element of the 1D array covered for each element (ties: the lower index).'''
    positions = np.nonzero(covered)[0]
    j = np.arange(len(covered))
    k = np.searchsorted(positions, j)
    before = positions[np.clip(k - 1, 0, len(positions) - 1)]
    after = positions[np.clip(k, 0, len(positions) - 1)]
    return np.where(np.abs(j - before) <= np.abs(after - j), before, after)


def fill_gaps(grid, under, dx, dz, verbose=True):
    '''
    Give each sample point not covered by any unit the unit of the nearest covered point: first along
    its row (gaps between shapes and shapes that do not reach the sides of the domain, as the drawn
    units are mostly horizontal), then along its column for the rows without any unit (strips at the
    top or the bottom of the domain, or between two shapes along the whole domain).
    Return grid and under filled and the largest distance filled [m].
    '''
    empty = grid < 0
    if not empty.any():
        return grid, under, 0.0
    if empty.all():
        raise SVGInputError("No sample point is covered by a unit: check that the shapes are inside the setup shape.")
    z_idx, x_idx = np.indices(grid.shape)
    rows = ~empty.all(axis=1)  # rows with at least one covered point
    src_x = x_idx.copy()
    for i in np.nonzero(rows)[0]:
        src_x[i] = nearest_covered(~empty[i])
    src_z = nearest_covered(rows)[z_idx]  # rows without units take the nearest row with units
    src_x = src_x[src_z, x_idx]
    distance = np.hypot((src_z - z_idx) * dz, (src_x - x_idx) * dx)
    if verbose:
        iz0, ix0 = np.argwhere(empty)[0]
        print(f"Warning: {empty.sum()} sample points not covered by any unit (first one at column {ix0}, depth "
              f"{(iz0 + 0.5) * dz / 1e3:.2f} km) took the unit of the nearest covered point, up to "
              f"{distance.max() / 1e3:.2f} km away. Check the gaps between shapes if this is too far.")
    return grid[src_z, src_x], under[src_z, src_x], distance.max()


def remove_slivers(grid, under, dz, Lz, min_thickness, unit_names, verbose=True):
    '''
    Remove the thin layers left by shapes that overlap a little (e.g. a crust drawn a bit below the base
    of the air leaves a thin layer of the mantle painted under it). A sliver is a layer of the stack (see
    build_stack) thinner than min_thickness [m] in every column, when its unit has a thicker layer in the
    stack: units thin everywhere (e.g. a decollement) and the thin ends of thick layers are kept.
    The points of the slivers take the unit right below them (above them at the bottom of the domain).
    Return grid and under without the slivers.
    '''
    if not min_thickness or min_thickness <= 0:
        return grid, under
    stack, unique, base = build_stack(grid, dz, Lz, unit_names)
    top = np.vstack([base[1:], np.zeros(base.shape[1])])
    thickness = (base - top).max(axis=1)  # largest thickness of each layer of the stack
    thick_units = {u for u, t in zip(stack, thickness) if t >= min_thickness}
    slivers = [k for k, (u, t) in enumerate(zip(stack, thickness)) if t < min_thickness and u in thick_units]
    if not slivers:
        return grid, under

    zdepth = (np.arange(grid.shape[0]) + 0.5) * dz
    slot = (zdepth[:, None, None] < base[None, 1:, :]).sum(axis=1)  # layer of the stack of each point
    sliver = np.isin(slot, slivers)
    removed = [(unit_names[stack[k]], np.nonzero(base[k] > top[k])[0], thickness[k]) for k in slivers]
    removed = [(name, cols[0], cols[-1], t) for name, cols, t in removed]

    grid, under = grid.copy(), under.copy()
    for ix in np.nonzero(sliver.any(axis=0))[0]:
        keep = np.nonzero(~sliver[:, ix])[0]
        iz = np.nonzero(sliver[:, ix])[0]
        k = np.searchsorted(keep, iz)
        src = keep[np.minimum(k, len(keep) - 1)]  # first point below (larger depth), or the last one above
        grid[iz, ix], under[iz, ix] = grid[src, ix], under[src, ix]
    if verbose:
        for name, c0, c1, t in removed:
            print(f"Warning: removed a sliver of '{name}' (columns {c0}-{c1}, up to {t / 1e3:.2f} km thick, "
                  f"thinner than min_thickness = {min_thickness / 1e3:g} km): it takes the unit below it.")
    return grid, under

def build_stack(grid, dz, Lz, unit_names):
    '''
    Return the stack of units (indices, bottom to top) and the base depth of each
    slot, an array (n_layers, Nx) in m (positive downwards).

    Depths are positive downwards on purpose: the sample grid comes from the SVG, whose y axis grows
    downwards (row 0 of grid is the top of the domain), and the whole reader works in depth. The
    conversion to the convention of Mandyoc (z positive upwards, from -Lz to 0, row 0 at the base)
    happens once, in layers_to_builder (z = -depth), which creates the interfaces of MandyocBuilder.
    '''
    Nzs, Nx = grid.shape
    if (grid < 0).any():
        iz, ix = np.argwhere(grid < 0)[0]
        n = (grid < 0).sum()
        raise SVGInputError(f"{n} sample points are not covered by any unit "
                 f"(first one at column {ix}, depth {(iz + 0.5) * dz / 1e3:.2f} km). "
                 "Check gaps between shapes or draw a background rectangle.")

    # runs of each column, from bottom to top: (units, depth of the top of each run)
    columns = []
    for ix in range(Nx):
        col = grid[::-1, ix]
        change = np.nonzero(np.diff(col))[0]  # last index of each run (except the top one)
        units = tuple(int(u) for u in col[np.r_[change, Nzs - 1]])
        tops = [Lz - (c + 1) * dz for c in change] + [0.0]
        columns.append((units, tops))

    unique = sorted({c[0] for c in columns}, key=len, reverse=True)
    stack = list(unique[0])
    for seq in unique[1:]:
        stack = scs(stack, seq)

    slots_of = {seq: embed(seq, stack) for seq in unique}
    base = np.empty((len(stack), Nx))
    base[0] = Lz
    for ix, (units, tops) in enumerate(columns):
        slots = slots_of[units]
        r = -1
        for k in range(1, len(stack)):
            # base of slot k = top of the highest run assigned to a slot below k
            while r + 1 < len(slots) and slots[r + 1] < k:
                r += 1
            base[k, ix] = tops[r] if r >= 0 else Lz

    # Slots absent in a column lie at the top of the run below them. When that run is of the same unit as a
    # slot above them (e.g. the lower crust split by a seed), the run is shared by both slots and the absent
    # slots go to the depth of their top in the nearest column where they are present: the interfaces of a
    # body continue flat around it, instead of rising to the top of the unit around it.
    drawn = base.copy()
    present = drawn > np.vstack([drawn[1:], np.zeros((1, Nx))])
    for ix, (units, tops) in enumerate(columns):
        slots = slots_of[units]
        for r, a in enumerate(slots):
            above = slots[r + 1] if r + 1 < len(slots) else len(stack)
            b = next((k for k in range(a + 2, above) if stack[k] == stack[a]), None)
            if b is None:
                continue
            cols = np.nonzero(present[a + 1:b].any(axis=0))[0]
            if not cols.size:
                continue
            near = cols[np.argmin(np.abs(cols - ix))]
            base[a + 1:b + 1, ix] = np.clip(drawn[b, near], tops[r], drawn[a, ix])  # inside the run
    return stack, unique, base

def rescale_thickness(stack, base, grid, under, unit_names, units, Lz, dz):
    '''
    Scale the drawing so that the units with h get that reference thickness. The drawing is rescaled, not
    replaced: a layer with h is stretched by h / (its most common drawn thickness) in every column, so where it is
    drawn thinner or thicker it keeps that proportion (e.g. a crust thinned to half keeps half of h).

    Each point of the drawing gets a vertical scale factor and, in each column,
    the depth of the interfaces is integrated from the top:
      - layers (units that reach both sides of the domain) with h: h / reference, where the
        reference is their most common drawn thickness (where they are not
        thinned or thickened). Thinned regions keep their drawn proportion;
      - layers without h: 1 (drawn thickness);
      - bodies (units that do not reach both sides: seeds, basins, weak
        zones, ...) without h: the factor of the unit they are drawn over (host);
      - bodies with h: h / reference, and their host in the same column absorbs
        the difference, so the units below do not move.
    The bottom unit fills the rest of the domain, so it can not have h.

    Return the new base depths (n_layers, Nx) and a list of messages describing the scaling.
    '''
    Nzs, Nx = grid.shape
    n_units = len(unit_names)
    has_h = np.array([units[name]["thickness"] is not None for name in unit_names])
    present = np.array([(grid == u).any(axis=0) for u in range(n_units)])  # (n_units, Nx)
    is_body = present.any(axis=1) & ~(present[:, 0] & present[:, -1])  # layers reach both sides of the domain

    def host_of(u):
        '''Unit most covered by u'''
        hosts = under[(grid == u) & (under >= 0)]
        return int(np.bincount(hosts).argmax()) if hosts.size else -1

    # factor of the units with h
    factor = np.ones(n_units)
    report = []
    for u, name in enumerate(unit_names):
        if not has_h[u] or not present[u].any():
            continue
        if u == stack[0]:
            raise SVGInputError(f"Unit '{name}' is the bottom layer and fills the rest of the domain: remove h from its description.")
        thickness = (grid == u).sum(axis=0)[present[u]]  # drawn thickness in samples, where present
        values, counts = np.unique(thickness, return_counts=True)
        reference = values[np.argmax(counts)] * dz
        share = counts.max() / present[u].sum()
        if share < 0.2:
            print(f"Warning: unit '{name}' has no clear constant thickness (the most common one, "
                  f"{reference / 1e3:.1f} km, occurs in {100 * share:.0f}% of its columns). Check the scaling to h.")
        factor[u] = units[name]["thickness"] / reference
        kind = "body" if is_body[u] else "layer"
        report.append(f"{kind} '{name}': drawn {reference / 1e3:.1f} km (in {100 * share:.0f}% of its columns) "
                      f"-> h = {units[name]['thickness'] / 1e3:.1f} km (factor {factor[u]:.3f})")

    # effective factor of each unit: bodies without h follow their host
    def effective_factor(u, depth=0):
        if u < 0:
            return 1.0
        if has_h[u] or not is_body[u] or depth > n_units:
            return factor[u]
        return effective_factor(host_of(u), depth + 1)
    effective = np.array([effective_factor(u) for u in range(n_units)])

    def host_factor(mask):
        '''Factor of the units covered at the points of mask'''
        hosts = under[mask]
        return np.where(hosts >= 0, effective[np.maximum(hosts, 0)], 1.0)

    F = np.ones((Nzs, Nx))
    for u, name in enumerate(unit_names):
        mask = grid == u
        if is_body[u] and not has_h[u]:
            F[mask] = host_factor(mask)
            if mask.any() and (factor != 1).any():
                host = host_of(u)
                report.append(f"body '{name}': no h, follows the scaling of '{unit_names[host] if host >= 0 else '-'}'")
        else:
            F[mask] = factor[u]

    # bodies with h: their host absorbs, in each column, the difference to the host scaling
    for u, name in enumerate(unit_names):
        host = host_of(u)
        if not (is_body[u] and has_h[u]) or host < 0:
            continue
        mask = grid == u
        expected = np.zeros((Nzs, Nx))
        expected[mask] = host_factor(mask)
        diff = (factor[u] * mask.sum(axis=0) - expected.sum(axis=0)) * dz  # extra thickness of the body (m)
        host_mask = grid == host
        host_total = (F * host_mask).sum(axis=0) * dz
        cols = mask.any(axis=0) & (diff != 0)
        scale = np.ones(Nx)
        absorbed = cols & (host_total > 0)
        scale[absorbed] = (host_total[absorbed] - diff[absorbed]) / host_total[absorbed]
        n_short = np.sum(cols & ((host_total <= 0) | (scale < 0)))
        if n_short:
            print(f"Warning: in {n_short} column(s) '{unit_names[host]}' is too thin to absorb the thickness of "
                  f"'{name}' (h); the units below are pushed down there.")
        F[host_mask] *= np.broadcast_to(np.clip(scale, 0.0, None), (Nzs, Nx))[host_mask]

    # new depth of the boundaries between samples, integrated from the top
    depth_new = np.vstack([np.zeros((1, Nx)), np.cumsum(F * dz, axis=0)])
    new_base = np.take_along_axis(depth_new, np.round(base / dz).astype(int), axis=0)
    new_base[0] = Lz

    if new_base[1].max() > Lz:
        ix = int(np.argmax(new_base[1]))
        raise SVGInputError(f"The layers above the bottom unit are {new_base[1, ix] / 1e3:.1f} km thick at column {ix}, "
                 f"more than Lz = {Lz / 1e3:.1f} km. Reduce h of some units.")
    return new_base, report

## I think the `MandyocLayer` is equivalent to the interfaces within `MandyocBuilder`
## May this function should return a xr.Dataset to be loaded as the interfaces or
## the whole `MandyocBuilder.DTree.interfaces`. A delicated merge is necessary
def layers_from_svg(svg_file, Lx, Lz, Nx, dz=100.0, setup_id=None, min_thickness=None, verbose=True):
    '''
    Read the SVG and return the list of MandyocLayer (bottom to top) with their
    interfaces (base depth, m) set. Lx, Lz and dz in m.
    The region mapped to the setup is the shape labelled "setup" if there is one,
    the element setup_id if given, or else the page.

    Parameters
    ----------
    svg_file : str
        Path to the SVG file.
    Lx : float
        Horizontal length of the domain in meters.
    Lz : float
        Vertical length of the domain in meters.
    Nx : int
        Number of horizontal grid points.
    dz : float, optional
        Vertical grid spacing in meters. Default is 100.0.
    setup_id : str, optional
        ID of the SVG element to use as the domain. If None, the shape labelled "setup" is used if present, otherwise the entire page is used. Default is None.
    min_thickness : float, optional
        Thinner pieces of a unit that is thicker elsewhere are removed, in meters (see remove_slivers;
        0 keeps them). If None, min_thickness of the setup shape, or MIN_THICKNESS (1 km). Default is None.
    verbose : bool, optional
        If True, print information about the layers and scaling. Default is True.
    '''
    if setup_id is None or min_thickness is None:
        setup = read_setup(svg_file)
        if setup is not None:
            setup_id = setup_id or setup["id"]
            if min_thickness is None:
                min_thickness = setup.get("min_thickness")
    if min_thickness is None:
        min_thickness = MIN_THICKNESS
    shape_unit, units = read_units(svg_file)
    if not units:
        raise SVGInputError("No shape with description found in the SVG.")
    unit_names = list(units)

    x, zdepth, dz, grid, under, colors = rasterize(svg_file, shape_unit, unit_names, Lx, Lz, Nx, dz, setup_id)
    grid, under, _ = fill_gaps(grid, under, Lx / (Nx - 1), dz, verbose=verbose)
    grid, under = remove_slivers(grid, under, dz, Lz, min_thickness, unit_names, verbose=verbose)
    stack, unique, base = build_stack(grid, dz, Lz, unit_names)

    # consistency check: rebuild the unit field of the drawing from the interfaces (Mandyoc convention)
    rebuilt = np.array(stack)[(zdepth[:, None, None] < base[None, 1:, :].transpose(0, 2, 1)).sum(axis=2)]
    mismatch = (rebuilt != grid).sum()

    # scale the units with a reference thickness h
    base, scaling = rescale_thickness(stack, base, grid, under, unit_names, units, Lz, dz)

    # one MandyocLayer per unit, copied for each slot of the stack
    unit_layers = []
    for i, name in enumerate(unit_names):
        u = units[name]
        layer = MandyocLayer(name, u["rheology"], Nx=Nx, title=u["title"], color=colors.get(i, "#808080"),
                             reference_thickness=u["thickness"], **u["props"])
        unit_layers.append(layer)

    layers = []
    for k, iunit in enumerate(stack):
        layer = copy.copy(unit_layers[iunit])
        layer.interface = base[k].copy()
        layer.base_depth = base[k].max()
        layers.append(layer)

    if verbose:
        for line in scaling:
            print(f"Thickness scaling: {line}")
        print(f"Distinct column patterns: {len(unique)}")
        print("Layers (bottom to top):")
        for k, layer in enumerate(layers):
            print(f"  {k:2d}  {layer.layer_label:25s} {layer.rheology_name:12s} "
                  f"C={layer.params['C']:g} rho={layer.params['rho']:g} H={layer.params['H']:.4g}")
        print(f"n_interfaces = {len(layers) - 1}   # put this in param.txt")
        print(f"Sampling: dx = {x[1] - x[0]:.1f} m, dz = {dz:.1f} m; mismatched samples of the drawing after rebuild: {mismatch}")

    return layers


def layer_indices(layers, label):
    '''Return the indices (position in the stack) of the layers with the given label.'''
    return [i for i, layer in enumerate(layers) if layer.layer_label == label]

def slot_thickness(layers, i, ix):
    '''Thickness [m] of the i-th layer of the stack at column(s) ix.'''
    top = layers[i+1].interface[ix] if i < len(layers)-1 else 0.0
    return layers[i].interface[ix] - top

def layer_thickness(layers, label, ix):
    '''Total thickness [m] of the layers with the given label at column(s) ix.'''
    return sum(slot_thickness(layers, i, ix) for i in layer_indices(layers, label))


def layers_to_builder(layers, Lx, Lz, Nz=2, path=".", seed_label="weak_seed", verbose=False):
    '''
    MandyocBuilder (tapioca.inputs) with the interfaces and the materials of the layers (bottom to top,
    output of layers_from_svg), so that the generic builders of tapioca can be used with an SVG scenario.
    The interface ID i is the top of layers[i] and carries its material (layer.material()); the top layer
    (sticky air) is the ID -1. Nz only matters for the fields (create_materials_fields, ...).
    The layers labelled seed_label (weak zones) get weakening_seed = 1 (initial strain of 1.0), unless
    weakening_seed is given in their description.
    '''
    builder = MandyocBuilder(path, Nx=len(layers[0].interface), Nz=Nz, Lx=Lx, Lz=Lz, verbose=verbose)
    builder.create_interface(-1, 0.0, interface_name=f"{len(layers)-1:02d}_{layers[-1].layer_label}",
                             **layer_material(layers[-1], seed_label))
    for i, layer in enumerate(layers[:-1]):
        builder.create_interface(i, -1.0 * layers[i+1].interface, interface_name=f"{i:02d}_{layer.layer_label}",
                                 **layer_material(layer, seed_label))
    return builder


def layer_material(layer, seed_label="weak_seed"):
    '''
    Parameters of the interface of a layer (MandyocLayer.material). weakening_seed is a parameter of the
    description like the others; when it is not given, the layers labelled seed_label get 1 (initial strain).
    '''
    params = layer.material()
    if layer.layer_label == seed_label and "weakening_seed" not in layer.given:
        params["weakening_seed"] = 1.0
    return params


def write_interfaces(layers, Lx, Lz, fname="interfaces.txt"):
    '''
    Write the interfaces.txt of Mandyoc (parameters of the layers and the interfaces, bottom to top) with
    MandyocBuilder, to any file name: used by tapioca-svg-interfaces, which writes only the interfaces, without
    a scenario (MandyocBuilder.export_interfaces always writes <path>/interfaces.txt).
    '''
    with tempfile.TemporaryDirectory() as folder:  # export_interfaces always writes <path>/interfaces.txt
        layers_to_builder(layers, Lx, Lz, path=folder).export_interfaces(export="none")
        shutil.move(os.path.join(folder, "interfaces.txt"), os.path.abspath(fname))


def plot_layers(layers, Lx, Lz, fname):
    x = np.linspace(0, Lx, layers[0].Nx)
    fig, ax = plt.subplots(figsize=(14, 6))
    for k, layer in enumerate(layers):
        top = layers[k + 1].interface if k < len(layers) - 1 else np.zeros_like(x)
        ax.fill_between(x / 1e3, -layer.interface / 1e3, -top / 1e3, color=layer.color, lw=0)
    ax.set_xlim(0, Lx / 1e3)
    ax.set_ylim(-Lz / 1e3, 0)
    ax.set_xlabel("Distance [km]")
    ax.set_ylabel("z [km]")

    # one legend entry per unit, from top to bottom
    unique = {}
    for layer in layers[::-1]:
        unique.setdefault(layer.layer_label, layer)
    handles = [Patch(facecolor=layer.color, edgecolor="k",
                     label=f"{layer.title}\n{layer.params['C']:g} x {layer.rheology_name.replace('_', ' ')}\n"
                           fr"$\rho$ = {layer.params['rho']:g} kg/m³")
               for layer in unique.values()]
    fig.subplots_adjust(bottom=0.28)
    ax.legend(handles=handles, ncol=len(handles), loc="lower center", bbox_to_anchor=(0.5, -0.42),
              frameon=False, title="Layers properties", title_fontsize=12, fontsize=8)
    fig.savefig(fname, dpi=200)
    plt.close(fig)


###############################################################################
# Main
###############################################################################

def main():
    '''Command line interface: tapioca-svg-interfaces scenario.svg -o interfaces.txt --plot setup.png'''
    try:
        _main()
    except SVGInputError as err:
        sys.exit(f"Error: {err}")


def _main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("svg", help="SVG file with the scenario")
    parser.add_argument("--Lx", type=float, help="domain width [km] (default: from the setup shape)")
    parser.add_argument("--Lz", type=float, help="domain depth [km] (default: from the setup shape)")
    parser.add_argument("--dx", type=float, help="horizontal resolution [km] (default: from the setup shape)")
    parser.add_argument("--Nx", type=int, help="number of points in x (alternative to --dx)")
    parser.add_argument("--sampling", type=float, default=0.1, help="vertical sampling to locate the interfaces [km] (default 0.1)")
    parser.add_argument("--min-thickness", type=float, help="remove thinner pieces [km] of units drawn thicker elsewhere, "
                        "left by overlapping shapes; 0 keeps them (default: from the setup shape, or 1 km)")
    parser.add_argument("--setup-id", help="id of an element whose bounding box is the domain (default: the shape labelled 'setup', or the page)")
    parser.add_argument("-o", "--output", default="interfaces.txt")
    parser.add_argument("--plot", help="save a figure of the resulting layers (e.g. numerical_setup.png)")
    args = parser.parse_args()

    SVGChecks()._check_svg_versions(args.svg, *svg_versions(args.svg), SVG_VERSION, INKSCAPE_VERSION)
    setup = read_setup(args.svg) or {}
    Lx = args.Lx * 1.0e3 if args.Lx else setup.get("Lx")
    Lz = args.Lz * 1.0e3 if args.Lz else setup.get("Lz")
    if Lx is None or Lz is None:
        raise SVGInputError(f"No shape labelled '{SETUP_LABEL}' in the SVG: give --Lx and --Lz (and --dx or --Nx).")
    if args.dx:
        Nx = points_from_spacing(Lx, args.dx * 1.0e3, "dx")
    else:
        Nx = args.Nx or setup.get("Nx")
    if Nx is None:
        raise SVGInputError("Give the horizontal resolution with --dx or --Nx.")
    print(f"Domain: Lx = {Lx / 1e3:g} km, Lz = {Lz / 1e3:g} km, Nx = {Nx}")

    min_thickness = args.min_thickness * 1.0e3 if args.min_thickness is not None else None
    layers = layers_from_svg(args.svg, Lx, Lz, Nx, dz=args.sampling * 1.0e3, setup_id=args.setup_id, min_thickness=min_thickness)
    write_interfaces(layers, Lx, Lz, args.output)
    print(f"Wrote {args.output}")

    if args.plot:
        plot_layers(layers, Lx, Lz, args.plot)
        print(f"Wrote {args.plot}")


if __name__ == "__main__":
    main()
