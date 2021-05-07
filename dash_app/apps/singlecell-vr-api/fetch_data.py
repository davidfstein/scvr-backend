import urllib.request
import os.path

links = {"10xPBMC_Seurat.loom": "https://drive.google.com/uc?export=download&id=1drBbfwq2kEWgkTF9WI907OvkiZpz9l18",
 "Macosko2015_Scanpy.h5ad": "https://drive.google.com/uc?export=download&id=1nNe4sbdTpPL8wkXUK-EA00GQrKzL9JGK", 
 "Nestorowa16_Stream.pkl": "https://drive.google.com/uc?export=download&id=1ajELnwuzTMlLdhrrvluUN85DWKSluKPK",
#  "Pancrease_Velocity.h5ad": "https://drive.google.com/uc?export=download&confirm=aW09&id=1fLblZqZ9oQJN12yb-qb7fFMZFdvkdOKf",
 "Paul2015_Paga.h5ad": "https://drive.google.com/uc?export=download&id=1HBk8JfGEZMMkYjJWoCsAbdXDtncEZDcc"}

def download_data():
    print(os.getcwd())
    for name, url in links.items():
        if os.path.exists('./dash_app/apps/singlecell-vr-api/app_datasets/%s'%name):
            continue
        urllib.request.urlretrieve(url, filename='./dash_app/apps/singlecell-vr-api/app_datasets/%s'%name)
        print('Done')