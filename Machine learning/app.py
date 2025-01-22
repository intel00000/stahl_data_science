from dash import Dash, dcc, html, Input, Output, callback
import os, json
import plotly.graph_objects as go
import pandas as pd
import base64

random_seed = 42
data_folder = "data"

# load the dataframe from pickle files
select_properties_df = pd.read_pickle(
    os.path.join(data_folder, "Select_properties.pkl")
)
yields_df = pd.read_pickle(os.path.join(data_folder, "Yields.pkl"))
yield_data_df = pd.read_pickle(os.path.join(data_folder, "yield_data_df.pkl"))
select_properties_data_df = pd.read_pickle(
    os.path.join(data_folder, "select_properties_data_df.pkl")
)
select_properties_data_removed_highlycorr_df = pd.read_pickle(
    os.path.join(data_folder, "select_properties_data_removed_highlycorr_df.pkl")
)
# load the dict from the json file
with open(os.path.join(data_folder, "mol_image_paths.json"), "r") as f:
    mol_image_paths = json.load(f)
with open(os.path.join(data_folder, "mol_image_data.json"), "r") as f:
    mol_image_data = json.load(f)
with open(os.path.join(data_folder, "mol_image_paths_captioned.json"), "r") as f:
    mol_image_paths_captioned = json.load(f)
with open(os.path.join(data_folder, "mol_image_data_captioned.json"), "r") as f:
    mol_image_data_captioned = json.load(f)

# Build the heatmap figure
fig = go.Figure(
    data=go.Heatmap(
        z=yields_df[yield_data_df.columns].T.to_numpy(),
        x=yields_df["id"],
        y=yield_data_df.columns,
        colorscale="Reds",
        colorbar=dict(
            title="Yield",
            dtick=20,
            tickvals=[0, 20, 40, 60, 80, 100],
            ticktext=["0%", "20%", "40%", "60%", "80%", "100%"],
        ),
        text=yields_df[yield_data_df.columns].T.to_numpy(),
        texttemplate="%{text:.2f}",
        hoverinfo="none",
        hovertemplate=None,  # Turn off default hovertemplate
        showscale=True,
    )
)

# Dash app setup
app = Dash()

app.layout = html.Div(
    [
        html.H1(children="Hello Dash"),
        dcc.Graph(id="heatmap-graph", figure=fig, clear_on_unhover=True),
        dcc.Tooltip(id="heatmap-tooltip"),
    ]
)

@app.callback(
    Output("heatmap-tooltip", "show"),
    Output("heatmap-tooltip", "bbox"),
    Output("heatmap-tooltip", "children"),
    Input("heatmap-graph", "hoverData"),
)
def display_hover(hoverData):
    print(hoverData)
    if hoverData is None:
        return False, None, None

    pt = hoverData["points"][0]
    x = pt["x"]  # Compound ID
    y = pt["y"]  # Method
    yield_value = pt["z"]

    # Get the image and compound details
    img_path = mol_image_paths.get(x)
    if img_path:
        # img_src = f"data:image/png;base64,{base64.b64encode(mol_image_data[img_path]).decode()}"
        # give the path directly
        base_file = os.path.basename(img_path)
        img_src = f"assets/images/{base_file}"
    else:
        img_src = None

    print(img_src)

    bbox = pt["bbox"]
    children = [
        html.Div(
            [
                html.Img(src=img_src, style={"width": "100%"}) if img_src else None,
                html.H3(f"Compound ID: {x}", style={"color": "darkblue"}),
                html.P(f"Method: {y}"),
                html.P(f"Yield: {yield_value:.2f}%"),
            ],
            style={"width": "200px", "white-space": "normal"},
        )
    ]

    return True, bbox, children


if __name__ == "__main__":
    app.run(debug=True)
