import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def pasv_post(file):
    df = pd.read_csv(file, sep='\t')

    gapless = df[~df['signature'].str.contains('-')].copy()


    gapless['orf_length'] = gapless['name'].apply(lambda x: (abs(int(x.split('_')[2]) - int(x.split('_')[1])) + 1) // 3)


    def classify_span(s):
        s = s.lower()
        if 'both' in s:
            return 'both'
        elif 'start' in s:
            return 'start'
        elif 'end' in s:
            return 'end'
        elif 'neither' in s:
            return 'neither'

    gapless['span_class'] = gapless['spans'].apply(classify_span)


    tally = gapless.groupby(['signature', 'span_class']).size().reset_index(name='count')

    print("Tally of signature and span class combinations:")
    print(tally)

    return gapless, tally


def boxplots(df, tally, output_file, protein):
    span_classes = df['span_class'].unique()
    
    sig_counts = [
        tally[tally['span_class'] == span_class]['signature'].nunique()
        for span_class in span_classes
    ]
    width_ratios = [max(1, n) *0.8 for n in sig_counts]
    fig_width = sum(width_ratios)

    orf_min = df['orf_length'].min()
    orf_max = df['orf_length'].max()

    fig, axes = plt.subplots(nrows=1, ncols=len(span_classes), figsize=(fig_width, 8), sharey=True, gridspec_kw={'width_ratios': width_ratios})

    if len(span_classes) == 1:
        axes = [axes]
    
    for i, span_class in enumerate(span_classes):
        ax = axes[i] 
        sub = df[df['span_class'] == span_class]

        sub_tally = (
            tally[tally['span_class'] == span_class]
            .sort_values(by='count', ascending=False)
        )
        signature_sort = sub_tally['signature'].tolist()
        data = [sub[sub['signature'] == sig]['orf_length'] for sig in signature_sort]
        ax.boxplot(data, labels=signature_sort)
        

        # y_min, y_max = ax.get_ylim()
        # padding = (y_max - y_min) * 0.05
        # ax.set_ylim(y_min - padding, y_max + padding)

        ax.yaxis.set_tick_params(labelleft=True)
        ax.set_yticks(ax.get_yticks())

        for j, sig in enumerate(signature_sort):
            values = sub[sub['signature'] == sig]['orf_length']
            count = sub_tally[sub_tally['signature'] == sig]['count'].values
            mean_val = round(values.mean())
            ax.annotate(f"N={count[0] if len(count) > 0 else 0}",
                xy=(j + 1, 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -20), textcoords='offset points',
                ha='center', va='top', fontsize=10, color='blue')

            ax.annotate(f"Mean={mean_val}",
                        xy=(j + 1, values.max()),
                        xytext=(0, 5), textcoords='offset points',
                        ha='center', va='bottom', fontsize=9, color='green')
            
        ax.set_title(f"Span: {span_class}")
        ax.set_ylabel("ORF Length (aa)")


    fig.suptitle(f"Signature Frequency of {protein}", fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()




def main(file, output_file, protein):
    df, tally = pasv_post(file)
    boxplots(df, tally, output_file, protein)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    protein = sys.argv[3]
    main(input_file, output_file, protein)



# def pasv_post(file):
#     df = pd.read_csv(file, sep='\t')
#     gapless = df[~df['signature'].str.contains('-')].copy()

#     gapless['orf_length'] = gapless['name'].apply(lambda x: (abs(int(x.split('_')[2]) - int(x.split('_')[1])) + 1) // 3)

#     def classify_span(s):
#         s = s.lower()
#         if 'both' in s:
#             return 'both'
#         elif 'start' in s:
#             return 'start'
#         elif 'end' in s:
#             return 'end'
#         elif 'neither' in s:
#             return 'neither'

#     gapless['span_class'] = gapless['spans'].apply(classify_span)
#     tally = gapless.groupby(['signature', 'span_class']).size().reset_index(name='count')

#     print("Tally of signature and span class combinations:")
#     print(tally)

#     return gapless, tally

# def boxplots(df, tally, output_file, protein):
#     span_classes = df['span_class'].unique()
    
#     sig_counts = [
#         tally[tally['span_class'] == span_class]['signature'].nunique()
#         for span_class in span_classes
#     ]
#     width_ratios = [max(1, n) *0.8 for n in sig_counts]
#     fig_width = sum(width_ratios)

#     orf_min = df['orf_length'].min()
#     orf_max = df['orf_length'].max()

#     fig, axes = plt.subplots(nrows=1, ncols=len(span_classes), figsize=(fig_width, 8), sharey=True, gridspec_kw={'width_ratios': width_ratios})

#     if len(span_classes) == 1:
#         axes = [axes]
    
#     for i, span_class in enumerate(span_classes):
#         ax = axes[i] 
#         sub = df[df['span_class'] == span_class]

#         sub_tally = (
#             tally[tally['span_class'] == span_class]
#             .sort_values(by='count', ascending=False)
#         )
#         signature_sort = sub_tally['signature'].tolist()
#         data = [sub[sub['signature'] == sig]['orf_length'] for sig in signature_sort]
#         ax.boxplot(data, labels=signature_sort)

#         ax.yaxis.set_tick_params(labelleft=True)
#         ax.set_yticks(ax.get_yticks())

#         for j, sig in enumerate(signature_sort):
#             values = sub[sub['signature'] == sig]['orf_length']
#             count = sub_tally[sub_tally['signature'] == sig]['count'].values
#             mean_val = round(values.mean())
#             ax.annotate(f"N={count[0] if len(count) > 0 else 0}",
#                 xy=(j + 1, 0), xycoords=('data', 'axes fraction'),
#                 xytext=(0, -20), textcoords='offset points',
#                 ha='center', va='top', fontsize=10, color='blue')

#             ax.annotate(f"Mean={mean_val}",
#                         xy=(j + 1, values.max()),
#                         xytext=(0, 5), textcoords='offset points',
#                         ha='center', va='bottom', fontsize=9, color='green')
            
#         ax.set_title(f"Span: {span_class}")
#         ax.set_ylabel("ORF Length (aa)")

#     fig.suptitle(f"Signature Frequency of {protein}", fontsize=16)
#     plt.tight_layout()
#     plt.savefig(output_file)
#     plt.show()

# def main(input_file, output_file, protein):
#     print(f"Processing {protein} data from {input_file}")
#     df, tally = pasv_post(input_file)
#     boxplots(df, tally, output_file, protein)
#     print(f"Plot saved to {output_file}")

# if __name__ == "__main__":
#     if len(sys.argv) != 4:
#         print("Usage: python pasv_post.py <input_tsv> <output_png> <protein_name>")
#         sys.exit(1)
    
#     input_file = sys.argv[1]
#     output_file = sys.argv[2]
#     protein = sys.argv[3]
#     main(input_file, output_file, protein)