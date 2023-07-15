import numpy as np
import pandas as pd
from flask import Flask, render_template, request, jsonify

app = Flask(__name__, static_folder='static')

abund_file = '9606_abund.txt'
domain_file = '9606_gn_dom.txt'

abund_df = pd.read_csv(abund_file, delimiter='\t')
domain_df = pd.read_csv(domain_file, delimiter='\t')

merged_df = pd.merge(abund_df, domain_df, how='inner', on='Gn')


def calculate_unique_values(df):
    return {'answer': df['Gn'].nunique()}


def calculate_mean_std(df):
    if 'std' not in df.columns:
        df['std'] = np.nan

    df['Mean-copy-number'] = pd.to_numeric(df['Mean-copy-number'], errors='coerce')
    df['std'] = pd.to_numeric(df['std'], errors='coerce')

    df.loc[df['std'].isnull(), 'std'] = df['Mean-copy-number'].mean()
    df.dropna(subset=['Mean-copy-number'], inplace=True)

    result = df[['Taxid', 'Ensembl_protein', 'Gn', 'Mean-copy-number']].to_dict(orient='records')
    return {'results': result}


def calculate_percentiles(df):
    df['Mean-copy-number'] = pd.to_numeric(df['Mean-copy-number'], errors='coerce')
    df = df.dropna(subset=['Mean-copy-number'])

    percentiles = df.groupby('Gn')['Mean-copy-number'].rank(pct=True) * 100
    return {'percentiles': percentiles.to_dict()}


def calculate_max_domain(df):
    df['Mean-copy-number'] = pd.to_numeric(df['Mean-copy-number'], errors='coerce')
    df = df.dropna(subset=['Mean-copy-number'])

    domain_avg = df.groupby('Domain')['Mean-copy-number'].mean()
    max_domain = domain_avg.idxmax()
    return {'answer': max_domain}


def calculate_domain_mean_std(df):
    df['Mean-copy-number'] = pd.to_numeric(df['Mean-copy-number'], errors='coerce')

    if 'std' not in df.columns:
        df['std'] = np.nan

    df.loc[df['Mean-copy-number'] == 'none', 'Mean-copy-number'] = np.nan

    df.dropna(subset=['Mean-copy-number'], inplace=True)

    df.loc[df['std'] == 'none', 'std'] = np.nan

    df['std'] = pd.to_numeric(df['std'], errors='coerce')

    df['std'].fillna(df['Mean-copy-number'].mean(), inplace=True)

    domain_stats = df.groupby(['Gn', 'Domain'])['Mean-copy-number'].agg(mean='mean', std='first').reset_index()

    result = domain_stats.to_dict(orient='records')
    return {'result': result}


def format_json_response(data):
    response = jsonify(data)
    response.headers['Content-Type'] = 'application/json'
    return response


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/calculate', methods=['POST'])
def calculate():
    question = request.form['question']

    if question == 'A1':
        result = calculate_unique_values(abund_df)
    elif question == 'A2':
        result = calculate_mean_std(abund_df)
    elif question == 'A3':
        result = calculate_percentiles(abund_df)
    elif question == 'B1':
        result = calculate_max_domain(merged_df)
    elif question == 'B2':
        result = calculate_domain_mean_std(merged_df)
    else:
        return format_json_response({'error': 'Invalid question'})

    return format_json_response(result)


if __name__ == '__main__':
    app.run(debug=True)
