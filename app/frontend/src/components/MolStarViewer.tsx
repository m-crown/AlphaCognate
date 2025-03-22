import { useEffect, useRef } from "react";
import { createPluginUI } from "molstar/lib/mol-plugin-ui";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";
import "molstar/lib/mol-plugin-ui/skin/light.scss";

declare global {
  interface Window {
    molstar?: PluginUIContext;
  }
}

interface MolStarViewerProps {
  pdbUrl: string;
}

export function MolStarViewer({ pdbUrl }: MolStarViewerProps) {
  const parentRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    async function initMolStar() {
      if (!parentRef.current) return;

      // Dispose existing instance
      window.molstar?.dispose();

      // Ensure MolStar does not overflow its container
      const container = parentRef.current;
      container.style.position = "relative"; // Needed for proper rendering

      // Initialize MolStar
      window.molstar = await createPluginUI({
        target: container,
        render: renderReact18,
      });

      // Load the selected PDB file
      const data = await window.molstar.builders.data.download(
        { url: pdbUrl },
        { state: { isGhost: true } }
      );
      const trajectory = await window.molstar.builders.structure.parseTrajectory(data, "pdb");
      await window.molstar.builders.structure.hierarchy.applyPreset(trajectory, "default");
    }

    initMolStar();

    return () => {
      window.molstar?.dispose();
      window.molstar = undefined;
    };
  }, [pdbUrl]);

  return (
    <div
      ref={parentRef}
      style={{
        maxHeight: '80vh',
        maxWidth: '100%',
        overflow: 'auto',
        width: "100%",
        height: "80vh", // Adjust height as needed
        border: "1px solid black",
        position: "relative",
      }}
    />
  );
}
