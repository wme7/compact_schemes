# Python3 module: compactSchemes
"""Python module that reads the alpha-coefs from db files"""
import sqlite3 as db

def get_alphas(db_filename,nbs_id):
    """Return the list of alpha coeficients in table nbs(nbs_id)"""
    conn = db.connect(db_filename)
    conn.row_factory = db.Row
    cur = conn.cursor()

    cur.execute("SELECT * FROM nbs WHERE id=?", (nbs_id,))
    row = cur.fetchone()

    params = {}
    for k, v in zip(row.keys()[1:], row[1:]):
        params[k] = float(v)

    return params