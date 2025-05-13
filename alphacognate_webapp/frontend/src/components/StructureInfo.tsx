import { Card, Text, Anchor, Stack } from '@mantine/core';
import ToggleBestButton from './ToggleBestButton';

interface Structure {
  name: string;
  num_transplants: number;
  runtime: number;
  url: string;
}

interface Props {
  structure: Structure | null;
  best: boolean;
  setBest: (value: boolean) => void;
}

export function StructureInfo({ structure, best, setBest }: Props) {
  if (!structure) return null;

  return (
    <Card shadow="sm" padding="md" radius="md" withBorder>
      <Stack gap="xs">
        <Text fw={500}>Name: {structure.name}</Text>
        <Text>Runtime: {structure.runtime.toFixed(2)}s</Text>
        <Text>Transplants: {structure.num_transplants}</Text>
        <Anchor href={structure.url} target="_blank" rel="noopener noreferrer">
          Download CIF
        </Anchor>
        <ToggleBestButton onToggle={(newVal) => setBest(newVal)} defaultChecked={best} />
      </Stack>
    </Card>
  );
}
