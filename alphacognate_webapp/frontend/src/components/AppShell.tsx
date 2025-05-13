import { AppShell, Burger, Group, Button, Stack } from '@mantine/core';
import { useDisclosure } from '@mantine/hooks';
import classes from '../css/AppShell.module.css';
import { Outlet } from 'react-router-dom';
import { Link } from 'react-router-dom';
import { MantineLogo } from '@mantinex/mantine-logo';
import { Logo } from './Logo';

export function MobileNavbar() {
  const [opened, { toggle }] = useDisclosure();

  return (
    <AppShell
      header={{ height: 60 }}
      navbar={{ width: 300, breakpoint: 'sm', collapsed: { desktop: true, mobile: !opened } }}
      padding="md"
    >
      <AppShell.Header>
        <Group h="100%" px="md">
          <Burger opened={opened} onClick={toggle} hiddenFrom="sm" size="sm" />
          <Group justify="space-between" style={{ flex: 1 }}>
            <Logo/>
            <Group ml="xl" gap={2} visibleFrom="sm">
              <Button className={classes.control} component={Link} to="/">Home</Button>
              <Button className={classes.control} component={Link} to="/about">About</Button>
            </Group>
          </Group>
        </Group>
      </AppShell.Header>

      <AppShell.Navbar py="md" px={4}>
        <Stack gap="md">
            <Button className={classes.control} component={Link} to="/">Home</Button>
            <Button className={classes.control} component={Link} to="/about">About</Button>
        </Stack>
      </AppShell.Navbar>

      <AppShell.Main>
        <Outlet/>
      </AppShell.Main>
    </AppShell>
  );
}