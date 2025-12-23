#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
from logging import exception

import numpy as np
import pandas as pd
from plotnine import *
import argparse
import argcomplete
import textwrap
from mizani.breaks import extended_breaks


COLUMNS = ["kernel_time_ms","output_time_ms","total_time_ms","nb_nodes","nb_steps","nb_snapshot","nb_sismo","nb_snapshot_ppm",
           "nb_histo","nb_slice","snapshot_output_size","slice_output_size","histogram_output_size","ppm_output_size","sismo_output_size",
           "total_snapshot_time","total_slice_time","total_histogram_time","total_ppm_time","total_sismo_time"]

def check_column_exists(df):
    for col in COLUMNS:
        if col not in df.columns:
            print(f"WARNING: The column '{col}' is missing in the CSV file")
    return True

def get_unit_label(col):
    if col in ["snapshot_output_size","slice_output_size","histogram_output_size","ppm_output_size","sismo_output_size"]:
        return "MB"
    elif col in ["kernel_time_ms","output_time_ms","total_time_ms","total_snapshot_time","total_slice_time","total_histogram_time","total_ppm_time","total_sismo_time"]:
        return "milliseconds"
    else :
        return ""

def filter_df(df, args):


    # no filter

    if df.empty:
        print("No data left after filtering. Please check your filters.")
        exit(1)

    return df

def pretreatment(df, x, y):

    # converte to MB
    size_cols = ["snapshot_output_size","slice_output_size","histogram_output_size","ppm_output_size","sismo_output_size"]
    for col in size_cols:
        if col in df.columns:
            df[col] = df[col] / (1024 * 1024)

    return df

def constantCols(df):
    return df.columns[df.eq(df.iloc[0]).all()].tolist()

def notConstantCols(df):
    return df.columns[~df.eq(df.iloc[0]).all()].tolist()

def createLegend(df, args, is_heatmap=False):
    x = args.x
    y = args.y
    spe_legend = args.spe_legend
    row = args.row
    _col = args.col

    ### black list of columns not to include in legend

    black_list = [x, y] + ["kernel_time_ms","output_time_ms","total_time_ms","snapshot_output_size","slice_output_size","histogram_output_size","ppm_output_size","sismo_output_size","total_snapshot_time","total_slice_time","total_histogram_time","total_ppm_time","total_sismo_time"]



    if is_heatmap and args.heatmap:
        black_list.append(args.heatmap)

    legend_cols = notConstantCols(df)
    legend_cols = [col for col in legend_cols if col not in black_list]

    if spe_legend:
        for col in spe_legend:
            if col in legend_cols:
                legend_cols.remove(col)

    if row and row in legend_cols:
        legend_cols.remove(row)

    if _col and _col in legend_cols:
        legend_cols.remove(_col)

    df['legend'] = df[legend_cols].apply(lambda row: ', '.join(f"{col}={row[col]}" for col in legend_cols), axis=1)

def calculate_perfect_scaling_df(df, x_col, y_col, scaling_type):
    perfect_scaling_points = []
    for legend_group, group_df in df.groupby('legend'):
        if group_df.empty: continue
        baseline_x = group_df[x_col].min()
        baseline_y = group_df[group_df[x_col] == baseline_x][y_col].mean()
        if baseline_x == 0: continue

        scaling_df = pd.DataFrame({'x_vals': sorted(group_df[x_col].unique())})
        if scaling_type == 'inverse':
            scaling_df['y_vals'] = (baseline_y * baseline_x) / scaling_df['x_vals']
        else:  # 'direct'
            scaling_df['y_vals'] = (baseline_y / baseline_x) * scaling_df['x_vals']

        scaling_df.rename(columns={'x_vals': x_col, 'y_vals': y_col}, inplace=True)
        scaling_df['legend'] = legend_group
        perfect_scaling_points.append(scaling_df)

    if not perfect_scaling_points: return pd.DataFrame()
    return pd.concat(perfect_scaling_points, ignore_index=True)

def get_axis_label(col_name, custom_name=None):
    if custom_name: return custom_name
    unit = get_unit_label(col_name)
    return f"{col_name} ({unit})" if unit else col_name

def normal_plot(df, x, y, title, legend_pos, x_label=None, y_label=None):
    p = (ggplot(df, aes(x=x, y=y, color='legend', fill='legend'))
         + stat_summary(geom='ribbon', alpha=0.3, size=0.001)
         + stat_summary(fun_y=np.mean, geom='line', size=2)
         + stat_summary(fun_y=np.mean, geom='point', size=7, shape='o', color='white', stroke=0.5)
         + scale_color_discrete(limits=df['legend'].unique())
         + scale_fill_discrete(limits=df['legend'].unique())
         + labs(title=title, x=get_axis_label(x, x_label), y=get_axis_label(y, y_label))
         + theme_bw()
         + theme(
                panel_grid_major=element_line(linestyle='--', color='lightgray'),
                panel_grid_minor=element_line(linestyle='--', color='lightgray', size=0.25),
                legend_text=element_text(size=16), legend_title=element_text(size=18),
                axis_title_x=element_text(size=18), axis_title_y=element_text(size=18),
                axis_text_x=element_text(size=16), axis_text_y=element_text(size=16),
                plot_title=element_text(size=18), legend_position=legend_pos,
                strip_text_x=element_text(size=18), strip_text_y=element_text(size=18),
            )
         )
    return p

def special_legend_plot(df, x, y, title, spe_legend, legend_pos, x_label=None, y_label=None):
    p = None
    if len(spe_legend) == 0:
        p = ggplot(df, aes(x=x, y=y, color='legend', fill='legend'))
    elif len(spe_legend) == 1:
        p = ggplot(df, aes(x=x, y=y, color='legend', fill='legend', linetype=f'factor({spe_legend[0]})'))
    elif len(spe_legend) == 2:
        p = ggplot(df, aes(x=x, y=y, color='legend', fill='legend', linetype=f'factor({spe_legend[0]})', shape=f'factor({spe_legend[1]})'))
    else:
        raise ValueError("You can only specify up to 2 special legend columns.")

    p += guides(shape=guide_legend(override_aes={"color": "black"}))
    p = ( p
          + stat_summary(geom='ribbon', alpha=0.3, size=0.001)
          + stat_summary(fun_y=np.mean, geom='line', size=2)
          + stat_summary(fun_y=np.mean, geom='point', size=7, color='black' if len(spe_legend) >= 2 else 'white', stroke=0.5)
          + scale_color_discrete(limits=df['legend'].unique())
          + scale_fill_discrete(limits=df['legend'].unique())
          + labs(
                linetype=spe_legend[0] if len(spe_legend) > 0 else None,
                shape=spe_legend[1] if len(spe_legend) > 1 else None,
                title=title, x=get_axis_label(x, x_label), y=get_axis_label(y, y_label),
            )
          + theme_bw()
          + theme(
                panel_grid_major=element_line(linestyle='--', color='lightgray'),
                panel_grid_minor=element_line(linestyle='--', color='lightgray', size=0.25),
                legend_text=element_text(size=16), legend_title=element_text(size=18),
                axis_title_x=element_text(size=18), axis_title_y=element_text(size=18),
                axis_text_x=element_text(size=16), axis_text_y=element_text(size=16),
                plot_title=element_text(size=18), legend_position=legend_pos, legend_key_width=30,
                strip_text_x=element_text(size=18), strip_text_y=element_text(size=18),
            )
          )
    return p

def heatmap_plot(df, x, y, z, title, legend_pos, x_label=None, y_label=None, row=None, col=None):

    # We group by X, Y and Facets (row/col) and take the mean of Z
    group_cols = [x, y]
    if row: group_cols.append(row)
    if col: group_cols.append(col)

    # Remove duplicates if row/col are the same as x/y (edge case)
    group_cols = list(set(group_cols))

    # Calculate Mean
    df = df.groupby(group_cols, as_index=False)[z].mean()

    # ----------------------------------------------------

    # 1. Force X and Y to be Ordered Categoricals (Factors) to create a perfect grid
    if pd.api.types.is_numeric_dtype(df[x]):
        df[x] = pd.Categorical(df[x], categories=sorted(df[x].unique()), ordered=True)
    else:
        df[x] = pd.Categorical(df[x], categories=sorted(df[x].unique()), ordered=True)

    if pd.api.types.is_numeric_dtype(df[y]):
        df[y] = pd.Categorical(df[y], categories=sorted(df[y].unique()), ordered=True)
    else:
        df[y] = pd.Categorical(df[y], categories=sorted(df[y].unique()), ordered=True)

    # 2. Prepare Labels
    def format_z(val):
        try:
            if abs(val) >= 100 or val == 0:
                return f"{val:.0f}"
            elif abs(val) < 0.01:
                return f"{val:.2e}"
            else:
                return f"{val:.2f}"
        except:
            return str(val)

    df['z_label'] = df[z].apply(format_z)

    # 3. Calculate text color contrast
    z_min, z_max = df[z].min(), df[z].max()
    mid_point = z_min + (z_max - z_min) * 0.45
    df['text_color'] = df[z].apply(lambda val: 'white' if val < mid_point else 'black')

    p = (ggplot(df, aes(x=x, y=y, fill=z))
         + geom_tile(color="white", size=0.5)
         + geom_text(aes(label='z_label', color='text_color'), size=10, show_legend=False)
         + scale_color_manual(values={'white': 'white', 'black': 'black'})

         + scale_fill_cmap(name='viridis')
         + labs(title=title, x=get_axis_label(x, x_label), y=get_axis_label(y, y_label), fill=get_axis_label(z))
         + theme_bw()
         + theme(
                panel_grid=element_blank(),
                legend_text=element_text(size=12), legend_title=element_text(size=14),
                axis_title_x=element_text(size=18), axis_title_y=element_text(size=18),
                axis_text_x=element_text(size=14, rotation=45, hjust=1), axis_text_y=element_text(size=14),
                plot_title=element_text(size=18), legend_position=legend_pos,
                strip_text_x=element_text(size=18), strip_text_y=element_text(size=18),
            )
         )
    return p

def add_reference_lines(p, lines_args, ref_x, ref_y):
    if not lines_args: return p
    for line_def in lines_args:
        try:
            parts = line_def.split(':')
            if len(parts) != 3: continue
            axis, label, value = parts
            value = float(value)

            if axis.lower() == 'y':
                p += geom_hline(yintercept=value, linetype='dashed', color='red', alpha=0.7, size=1)
                p += geom_text(x=ref_x, y=value, label=label, color='red', size=16, fontweight='bold', ha='left', va='bottom', inherit_aes=False)
            elif axis.lower() == 'x':
                p += geom_vline(xintercept=value, linetype='dashed', color='blue', alpha=0.7, size=1)
                p += geom_text(x=value, y=ref_y, label=label, color='blue', size=16, fontweight='bold', ha='left', va='bottom', inherit_aes=False, angle=90)
        except ValueError:
            print(f"WARNING: Could not parse value in line definition '{line_def}'. Skipping.")
    return p

def main():
    parser = argparse.ArgumentParser(description='Plot data from a CSV file.')
    parser.add_argument('--x', type=str, required=True, help='Column name for x-axis')
    parser.add_argument('--y', type=str, required=True, help='Column name for y-axis')
    parser.add_argument('--input', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--output', type=str, required=True, help='Path to save the plot')
    parser.add_argument('--delete', nargs='+', help='List of columns to delete')
    parser.add_argument('--xscale', type=str, choices=['linear', 'log', 'log2'], default='linear', help='Scale for x-axis')
    parser.add_argument('--yscale', type=str, choices=['linear', 'log', 'log2'], default='linear', help='Scale for y-axis')
    parser.add_argument('--y0', action='store_true', help='Set y-axis to start at 0')
    parser.add_argument('--x0', action='store_true', help='Set x-axis to start at 0')
    parser.add_argument('--spe-legend', nargs='+', help="List of special columns, for own legend/style (max 2 columns, style applied in order, linestyle, marker)")
    parser.add_argument('--legend-position', type=str, choices=['bottom', 'top', 'left', 'right', 'none'], default='right', help='Position of the legend.')
    parser.add_argument('--facet', type=str, choices=['grid', 'wrap'], default='grid', help='Facet type for the plot (default: grid)')
    parser.add_argument('--row', type=str, help='Column name for separated plots on rows')
    parser.add_argument('--row-free', action='store_true', help='Use free scales for row facets')
    parser.add_argument('--col', type=str, help='Column name for separated plots on columns')
    parser.add_argument('--col-free', action='store_true', help='Use free scales for column facets')
    parser.add_argument('--perfect-scaling', action='store_true', help='Show a perfect scaling reference line based on the first data point of each curve.')
    parser.add_argument('--scaling-type', type=str, choices=['inverse', 'direct'], default='direct',
                        help="""Type of perfect scaling: 'inverse' (e.g., time vs cpus, y=k/x) or 
                                'direct' (e.g., speed vs cpus, y=k*x). Default is 'direct'.""")
    parser.add_argument('--title', type=str, help='Custom title for the plot')
    parser.add_argument('--line', nargs='+', help='Add reference lines (axis:label:value) ex : --line y:TargetTime1:50 y:TargetTime1:35 x:MaxSize:1000')
    parser.add_argument('--xname', type=str, help='Custom label for X axis')
    parser.add_argument('--yname', type=str, help='Custom label for Y axis')
    parser.add_argument('--heatmap', type=str, help='Enable heatmap mode (Z axis) ex: --x block_size --y sub_block_size --heatmap GFlop_per_s')

    ####################### filter options ###########################


    # no filter

    ############################################################################

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    if args.spe_legend and len(args.spe_legend) > 2:
        raise ValueError("You can only specify up to 2 special legend columns.")

    try :
        df = pd.read_csv(args.input, delimiter=',')
    except FileNotFoundError:
        print(f"ERROR: The file '{args.input}' was not found.")
        exit(1)

    check_column_exists(df)
    df = pretreatment(df, args.x, args.y)
    df = filter_df(df, args)
    if args.delete: df = df.drop(columns=args.delete)

    # Check if x, y (and heatmap) columns exist
    for col in [args.x, args.y] + ([args.heatmap] if args.heatmap else []):
        if col not in df.columns:
            print(f"ERROR: The specified column '{col}' does not exist in the data.\nNormally available columns are: {', '.join(df.columns)}")
            exit(1)

    # TITLE LOGIC
    if args.title:
        wrapped_title = textwrap.fill(args.title, width=50)
    else:
        title_parts = []
        for col in constantCols(df):
            if len(df[col].unique()) == 1: title_parts.append(f"{col}={df[col].unique()[0]}")
        title_parts.append(f"nb_xp={len(df)}")
        title = ', '.join(title_parts)
        wrapped_title = textwrap.fill(title, width=50)

    createLegend(df, args, is_heatmap=bool(args.heatmap))

    # PLOT GENERATION
    if args.heatmap:
        print(f"Generating Heatmap: X={args.x}, Y={args.y}, Z={args.heatmap}")
        if not pd.api.types.is_numeric_dtype(df[args.heatmap]):
            print(f"WARNING: Heatmap column '{args.heatmap}' is not numeric.")
        # Pass ROW and COL to allow aggregation
        p = heatmap_plot(df, args.x, args.y, args.heatmap, wrapped_title, args.legend_position,
                         x_label=args.xname, y_label=args.yname, row=args.row, col=args.col)
    else:
        if args.spe_legend is None:
            p = normal_plot(df, args.x, args.y, wrapped_title, args.legend_position, x_label=args.xname, y_label=args.yname)
        else:
            p = special_legend_plot(df, args.x, args.y, wrapped_title, args.spe_legend, args.legend_position, x_label=args.xname, y_label=args.yname)

    # --- AXIS SCALING AND LIMITS (ONLY FOR NORMAL PLOTS) ---
    if not args.heatmap:
        try:
            x_min = 0 if args.x0 else df[args.x].min()
            x_max = df[args.x].max()
            y_min = 0 if args.y0 else df[args.y].min()
            y_max = df[args.y].max()
        except:
            x_min, x_max, y_min, y_max = 0, 1, 0, 1

        if args.line:
            for line_def in args.line:
                try:
                    parts = line_def.split(':')
                    if len(parts) == 3:
                        axis, _, val = parts
                        val = float(val)
                        if axis.lower() == 'y':
                            if val > y_max: y_max = val
                            if val < y_min: y_min = val
                        elif axis.lower() == 'x':
                            if val > x_max: x_max = val
                            if val < x_min: x_min = val
                except: pass

        exact_x_breaks = True
        exact_y_breaks = True
        x_breaks_values = sorted(df[args.x].unique().tolist()) if pd.api.types.is_numeric_dtype(df[args.x]) else []
        y_breaks_values = sorted(df[args.y].unique().tolist()) if pd.api.types.is_numeric_dtype(df[args.y]) else []

        if len(x_breaks_values) > 20: exact_x_breaks = False
        if len(y_breaks_values) > 20: exact_y_breaks = False

        if exact_x_breaks and pd.api.types.is_numeric_dtype(df[args.x]):
            filtered = [x_breaks_values[0]]
            last = x_breaks_values[0]
            denom = x_max if x_max!=0 else 1
            for val in x_breaks_values[1:]:
                if abs(val - last)/denom > 0.03:
                    filtered.append(val)
                    last = val
            if abs(filtered[-1]-x_max)/denom > 0.03 and x_max not in filtered: filtered.append(x_max)
            x_breaks_values = filtered

        if exact_y_breaks and pd.api.types.is_numeric_dtype(df[args.y]):
            filtered = [y_breaks_values[0]]
            last = y_breaks_values[0]
            denom = y_max if y_max!=0 else 1
            for val in y_breaks_values[1:]:
                if abs(val - last)/denom > 0.03:
                    filtered.append(val)
                    last = val
            if abs(filtered[-1]-y_max)/denom > 0.03 and y_max not in filtered: filtered.append(y_max)
            y_breaks_values = filtered

        x_limits = (x_min, x_max) if pd.api.types.is_numeric_dtype(df[args.x]) else None
        y_limits = (y_min, y_max) if pd.api.types.is_numeric_dtype(df[args.y]) else None
        if args.col or args.row: x_limits, y_limits = None, None

        trans_x = "log2" if args.xscale == "log2" else ("log10" if args.xscale == "log" else None)
        trans_y = "log2" if args.yscale == "log2" else ("log10" if args.yscale == "log" else None)

        if args.x0 and trans_x: print("WARNING: --x0 ignored because x-axis uses a logarithmic scale.")
        if args.y0 and trans_y: print("WARNING: --y0 ignored because y-axis uses a logarithmic scale.")

        if pd.api.types.is_numeric_dtype(df[args.x]):
            if trans_x: p += scale_x_continuous(trans=trans_x)
            elif exact_x_breaks: p += scale_x_continuous(breaks=x_breaks_values, limits=x_limits)
            else: p += scale_x_continuous(breaks=extended_breaks(n=10), limits=x_limits)

        if pd.api.types.is_numeric_dtype(df[args.y]):
            if trans_y: p += scale_y_continuous(trans=trans_y)
            elif exact_y_breaks: p += scale_y_continuous(breaks=y_breaks_values, limits=y_limits)
            else: p += scale_y_continuous(breaks=extended_breaks(n=10), limits=y_limits)

        if not trans_x and not trans_y:
            if args.x0 and args.y0 and (args.row or args.col): p += expand_limits(x=0, y=0)
            elif args.x0: p += expand_limits(x=0)
            elif args.y0: p += expand_limits(y=0)

    # FACETS
    if args.row or args.col:
        scales = 'fixed'
        if args.row_free and args.col_free: scales = 'free'
        elif args.row_free: scales = 'free_y'
        elif args.col_free: scales = 'free_x'
        facets = (args.row, args.col) if (args.row and args.col) else (args.row if args.row else args.col)
        nrow = df[args.row].nunique() if args.row else 1
        ncol = df[args.col].nunique() if args.col else 1
        if args.facet == 'grid': p += facet_grid(rows=args.row, cols=args.col, scales=scales, labeller="label_both")
        else: p += facet_wrap(facets=facets, scales=scales, nrow=nrow, ncol=ncol, labeller="label_both")

    if args.perfect_scaling and not args.heatmap:
        print(f"Calculating perfect scaling line (type: {args.scaling_type})...")
        perfect_df = calculate_perfect_scaling_df(df, args.x, args.y, args.scaling_type)
        if not perfect_df.empty:
            p += geom_line(data=perfect_df, mapping=aes(x=args.x, y=args.y, color='legend'), linetype='dashed', size=1, inherit_aes=False)

    # ADD REF LINES
    if args.line:
        safe_x = df[args.x].min() if pd.api.types.is_numeric_dtype(df[args.x]) else 0
        safe_y = df[args.y].min() if pd.api.types.is_numeric_dtype(df[args.y]) else 0

        if not args.heatmap:
            try:
                x_min = 0 if args.x0 else df[args.x].min()
                x_max = df[args.x].max()
                y_min = 0 if args.y0 else df[args.y].min()
                y_max = df[args.y].max()
                if 'trans_x' in locals() and trans_x and safe_x <= 0:
                    pos_x = df[args.x][pd.to_numeric(df[args.x], errors='coerce') > 0]
                    safe_x = pos_x.min() if not pos_x.empty else 1e-9
                else: safe_x = x_min
                if 'trans_y' in locals() and trans_y and safe_y <= 0:
                    pos_y = df[args.y][pd.to_numeric(df[args.y], errors='coerce') > 0]
                    safe_y = pos_y.min() if not pos_y.empty else 1e-9
                else: safe_y = y_min
            except: pass

        p = add_reference_lines(p, args.line, safe_x, safe_y)

    p.save(args.output, width=16, height=9, units='in', dpi=300, verbose=False)
    print(f"Plot saved to {args.output}")
    print(f"try: open {args.output}")

if __name__ == "__main__":
    main()