import { useParams } from "react-router-dom";
import StructureTable from "../components/StructureTable";
import AutoComplete from "../components/SearchBox";
import { MolStarViewer } from "../components/MolStarViewer";
import { Grid } from "@mantine/core";

function StructurePage() {
  const { id } = useParams(); // Extract the ID from the URL

  // Here, make API call to get the necessary info for Table view and structure loading.
  // In the future, on click of table row, would be good to load a comparison of ligands.

  return (
        <Grid>
          <Grid.Col span={6}>
            <AutoComplete />
            <StructureTable />
          </Grid.Col>
          <Grid.Col span={6}>
            <MolStarViewer pdbUrl={id || ""} />
          </Grid.Col>
        </Grid>
  );
}

export default StructurePage;