#!/usr/bin/env python3
import pandas as pd
import plotly.graph_objects as go
import argparse
import glob
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Analyze Bracken files and generate interactive HTML.")
    parser.add_argument('files', nargs='+', help="Bracken input files (supports wildcards)")
    parser.add_argument('--output', '-o', default='bracken_analysis.html', help="Output HTML file")
    args = parser.parse_args()

    # Expand wildcards (just in case shell didn't)
    file_list = []
    for f in args.files:
        expanded = glob.glob(f)
        if expanded:
            file_list.extend(expanded)
        else:
            file_list.append(f)

    file_list = sorted(list(set(file_list)))

    if not file_list:
        print("No files found.")
        return

    dfs = []
    for f in file_list:
        try:
            # Try reading with tab separator first, handling comments
            try:
                df = pd.read_csv(f, sep='\t', comment='#')
            except Exception:
                # If that fails, try whitespace separator
                try:
                    df = pd.read_csv(f, sep=r'\s+', comment='#')
                except Exception:
                     # If both fail, raise the original error or a generic one
                     raise ValueError("Could not parse file with tab or whitespace separator.")

            # Check if we got a reasonable dataframe (more than 1 column)
            if df.shape[1] < 2:
                 # Retry with engine='python' which is more robust for some separators
                 df = pd.read_csv(f, sep=None, engine='python', comment='#')

            df.columns = df.columns.str.strip()

            required_cols = ['name', 'kraken_assigned_reads', 'added_reads', 'new_est_reads', 'fraction_total_reads']
            missing = [col for col in required_cols if col not in df.columns]
            if missing:
                 print(f"Skipping {f}: missing columns {missing}. Found: {list(df.columns)}", file=sys.stderr)
                 continue

            df['Source'] = os.path.basename(f)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}", file=sys.stderr)

    if not dfs:
        print("No valid data loaded.")
        return

    combined_df = pd.concat(dfs, ignore_index=True)

    metrics = ['kraken_assigned_reads', 'added_reads', 'new_est_reads', 'fraction_total_reads']
    metrics = [m for m in metrics if m in combined_df.columns]

    if not metrics:
        print("No numeric metrics found to plot.")
        return

    initial_metric = metrics[0]

    def get_sorted_species(metric):
        total_per_species = combined_df.groupby('name')[metric].sum().sort_values(ascending=False)
        return total_per_species.index.tolist()

    fig = go.Figure()

    sources = sorted(combined_df['Source'].unique())
    sorted_species_init = get_sorted_species(initial_metric)

    for source in sources:
        source_df = combined_df[combined_df['Source'] == source]
        source_df_indexed = source_df.set_index('name').reindex(sorted_species_init).fillna(0)

        fig.add_trace(go.Bar(
            x=source_df_indexed.index,
            y=source_df_indexed[initial_metric],
            name=source
        ))

    buttons = []
    for metric in metrics:
        sorted_species = get_sorted_species(metric)

        new_x = []
        new_y = []

        for source in sources:
            source_df = combined_df[combined_df['Source'] == source]
            source_df_indexed = source_df.set_index('name').reindex(sorted_species).fillna(0)

            new_x.append(source_df_indexed.index.tolist())
            new_y.append(source_df_indexed[metric].tolist())

        button = dict(
            label=metric,
            method="update",
            args=[{"x": new_x, "y": new_y},
                  {"yaxis.title.text": metric, "title": f"Bracken Analysis: {metric}"}]
        )
        buttons.append(button)

    fig.update_layout(
        updatemenus=[dict(
            type="buttons",
            direction="down",
            buttons=buttons,
            showactive=True,
            x=1.02,
            xanchor="left",
            y=1,
            yanchor="top"
        )],
        title=f"Bracken Analysis: {initial_metric}",
        xaxis_title="Species",
        yaxis_title=initial_metric,
        barmode='group',
        margin=dict(r=150)
    )

    fig.write_html(args.output)
    print(f"Created {args.output}")

if __name__ == "__main__":
    main()
