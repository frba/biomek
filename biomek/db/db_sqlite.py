import os
import sqlite3
from sqlite3 import Error

db_file = os.path.join('/home/flavia/Documents/Concordia/projects/Web_tools/', 'db.sqlite3')


def execute():
    conn = connect_db()
    if conn is not None:
        new_sample = ('MSH_Sdi1', '123456', '4766', '')
        add_sample(conn, 'sample', new_sample)
        # remove_sample(conn, 'sample', 'MSH_Sdi1')

        # rows = get_sample(conn, 'sample', 'MSH_Sdi1')
        rows = get_table(conn, 'sample')
        for row in rows:
            print(row)


def connect_db():
    try:
        conn = sqlite3.connect(db_file)
        # print("connect successful!")
        return conn
    except Error as e:
        print(e)
    return None


def get_sample(conn, table, sample_name):
    sql_show_samples = "SELECT * FROM "+str(table)+" WHERE name = '" + str(sample_name) + "'"
    cur = conn.cursor()
    cur.execute(sql_show_samples)
    rows = cur.fetchall()
    return rows


def get_table(conn, table):
    sql_show_samples = "SELECT * FROM "+str(table)+""
    cur = conn.cursor()
    cur.execute(sql_show_samples)
    rows = cur.fetchall()
    return rows


def add_sample(conn, table, new_sample):
    sql_insert_sample = " INSERT INTO "+str(table)+" ('name', 'type', 'length', 'sequence') VALUES (?,?,?,?)"
    cur = conn.cursor()
    cur.execute(sql_insert_sample, new_sample)
    conn.commit()


def remove_sample(conn, table, sample_name,):
    sql_remove_sample = " DELETE from "+str(table)+" WHERE name = '" + str(sample_name) + "'"
    cur = conn.cursor()
    cur.execute(sql_remove_sample)
    conn.commit()
