import pymysql, sys
# from pydent import AqSession
import pydent as pd


# def connect_db():
#     try:
#         connection = pymysql.connect(host='172.19.0.2',
#                                  user='root',
#                                  password='aSecretAquarium',
#                                  db='production',
#                                  charset='utf8mb4',
#                                  cursorclass=pymysql.cursors.DictCursor)
#     except:
#         try:
#             connection = pymysql.connect(host='172.19.0.3',
#                                          user='root',
#                                          password='aSecretAquarium',
#                                          db='production',
#                                          charset='utf8mb4',
#                                          cursorclass=pymysql.cursors.DictCursor)
#         except:
#             print('Not possible to connect to the database')
#             sys.exit()
#
#     print("connect successful!!")
#     return connection


def connect_db():
    try:
        connection = pd.AqSession('flavia', 'genomefoundry', 'http://flavia-optiplex-3060:8080/')
    except:
        try:
            connection = pd.AqSession('root', 'aSecretAquarium', '172.19.0.2')
        except:
            print('Not possible to connect to the database')
            sys.exit()

    sample = connection.Sample.find('ConR1')
    print(sample)
    print("connect successful!!")
    return connection


def get_sample(connection, name):
    with connection.cursor() as cursor:
        sql = "SELECT samples.id, samples.name, sample_type_id, project, user_id, description " \
              "FROM production.samples WHERE samples.name = "+"'"+ str(name)+"'"
        cursor.execute(sql)
        if cursor.rowcount > 0:
            for row in cursor:
                print(row)
            return cursor
        else:
            print('Sample: ' + str(name) + ' not found!')
            return None


def get_samples_from_project(connection, project):
    with connection.cursor() as cursor:
        sql = "SELECT samples.id, samples.name, sample_type_id, project, user_id, description " \
              "FROM production.samples WHERE project = "+"'"+ str(project)+"'"
        cursor.execute(sql)
        for row in cursor:
            print(row)
    return cursor


def get_sample_in_plate(connection, name):
    with connection.cursor() as cursor:
        sql = "SELECT samples.id, samples.name, sample_types.name, " \
              "items.sample_id, items.id, " \
              "part_associations.collection_id, part_associations.row, part_associations.column, " \
              "data_associations.object " \
              "FROM production.samples, production.sample_types, production.items, " \
              "production.part_associations, production.data_associations " \
              "WHERE samples.name = " + "'" + str(name) + "'" + "AND sample_type_id = sample_types.id " \
                                                                "AND samples.id = items.sample_id " \
                                                                "AND items.quantity > 0 " \
                                                                "AND items.id = part_associations.part_id " \
                                                                "AND items.id = data_associations.parent_id"
        cursor.execute(sql)
        for row in cursor:
            print(row)
    return cursor


def get_sample_info(connection, name):
    with connection.cursor() as cursor:
        sql = "SELECT samples.id, samples.name, sample_types.name, items.sample_id, items.id, field_values.name, field_values.value " \
              "FROM production.samples, production.sample_types, production.items, production.field_values " \
              "WHERE samples.name = " + "'" + str(name) + "'" + "AND sample_type_id = sample_types.id " \
                                                                "AND samples.id = items.sample_id " \
                                                                "AND items.quantity > 0 " \
                                                                "AND samples.id = field_values.parent_id "
        cursor.execute(sql)
        for row in cursor:
            print(row)
    return cursor


def get_info(connection):
    with connection.cursor() as cursor:
        sql = "select user,host from mysql.user;"
        cursor.execute(sql)
        for row in cursor:
            print(row)
    return cursor


def get_data(connection):
    with connection.cursor() as cursor:
        # SQL
        sql = "SELECT * FROM production.object_types WHERE name = '384' "

        # Execute query.
        cursor.execute(sql)

        # print("cursor.description: ", cursor.description)
        #
        # print()

        for row in cursor:
            print(row)

    return cursor


def insert(connection):
    with connection.cursor() as cursor:
        sql = "INSERT INTO production.users VALUES (5,'Audrey','audrey','2013-06-16 17:26:54','2017-10-19 04:59:18','$2a$10$HxgxLX5/ITcYpII1InAL1.jUYAiHk/rMftHniPJVvauy43VDoo8yW','TYmoWfyV42AL7dSoYcgmug',1,'VHzz9IW3xnNx8O3cA_P0rKsUWmTVH_Qz9mHKqgE-hNI');"
        cursor.execute(sql)
        connection.commit()
    return cursor


def update(connection):
    with connection.cursor() as cursor:
        sql = "UPDATE production.users SET name = 'SMITA' WHERE id = 4"
        cursor.execute(sql)
        connection.commit()
    return cursor





